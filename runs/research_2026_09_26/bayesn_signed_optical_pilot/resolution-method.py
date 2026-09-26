def _setup_band_weights(self):
    """
    Sets up the interpolation for the band weights used for photometry as well as calculating the zero points for
    each band. This code is partly based off ParSNiP from Boone+21
    """
    # Build the model in log wavelength
    self.min_wave = self.l_knots[0]
    self.max_wave = self.l_knots[-1]
    self.spectrum_bins = self.audit_spectrum_bins
    self.band_oversampling = 51
    self.max_redshift = 4

    model_log_wave = np.linspace(np.log10(self.min_wave),
                                 np.log10(self.max_wave),
                                 self.spectrum_bins)

    model_spacing = model_log_wave[1] - model_log_wave[0]

    band_spacing = model_spacing / self.band_oversampling
    band_max_log_wave = (
            np.log10(self.max_wave * (1 + self.max_redshift))
            + band_spacing
    )

    # Oversampling must be odd.
    assert self.band_oversampling % 2 == 1
    pad = (self.band_oversampling - 1) // 2
    band_log_wave = np.arange(np.log10(self.min_wave),
                              band_max_log_wave, band_spacing)
    band_wave = 10 ** band_log_wave

    # Load in-built filter yaml first
    with open(os.path.join(self.__root_dir__, 'bayesn-filters', 'filters.yaml'), 'r') as file:
        filter_dict = yaml.load(file)

    # Prepend root locations for in-built filters
    for key, val in filter_dict['standards'].items():
        filter_dict['standards'][key]['path'] = os.path.join(self.__root_dir__, 'bayesn-filters', val['path'])

    for key, val in filter_dict['filters'].items():
        filter_dict['filters'][key]['path'] = os.path.join(self.__root_dir__, 'bayesn-filters', val['path'])

    # Add custom filters, if specified
    if self.filter_yaml is not None:
        if not os.path.exists(self.filter_yaml):
            raise FileNotFoundError(f'Specified filter yaml {self.filter_yaml} does not exist')
        with open(self.filter_yaml, 'r') as file:
            custom_filter_dict = yaml.load(file)
        # Add custom standards if specified---------------------
        if 'standards' in custom_filter_dict.keys():
            if 'standards_root' in custom_filter_dict.keys():
                standards_root = custom_filter_dict['standards_root']
            else:
                standards_root = ''
            for key, val in custom_filter_dict['standards'].items():
                path = os.path.join(standards_root, val['path'])
                # Fill environment variables if used e.g. $SNDATA_ROOT
                split_path = os.path.normpath(path).split(os.path.sep)
                root = split_path[0]
                if root[:1] == '$':
                    env = os.getenv(root[1:])
                    if env is None:
                        raise FileNotFoundError(f'The environment variable {root} was not found')
                    path = os.path.join(env, *split_path[1:])
                elif not os.path.isabs(path):  # If relative path, prepend yaml location
                    path = os.path.join(os.path.split(os.path.abspath(self.filter_yaml))[0], path)
                custom_filter_dict['standards'][key]['path'] = path
                # Add custom standard and overwrite existing one of same name if present
                filter_dict['standards'][key] = custom_filter_dict['standards'][key]
        # Add custom filters
        if 'filters_root' in custom_filter_dict.keys():
            filters_root = custom_filter_dict['filters_root']
        else:
            filters_root = ''
        for key, val in custom_filter_dict['filters'].items():
            path = os.path.join(filters_root, val['path'])
            # Fill environment variables if used e.g. $SNDATA_ROOT
            split_path = os.path.normpath(path).split(os.path.sep)
            root = split_path[0]
            if root[:1] == '$':
                env = os.getenv(root[1:])
                if env is None:
                    raise FileNotFoundError(f'The environment variable {root} was not found')
                path = os.path.join(env, *split_path[1:])
            elif not os.path.isabs(path):  # If relative path, prepend yaml location
                path = os.path.join(os.path.split(os.path.abspath(self.filter_yaml))[0], path)
            custom_filter_dict['filters'][key]['path'] = path
            # Add custom filter and overwrite existing one of same name if present
            filter_dict['filters'][key] = custom_filter_dict['filters'][key]

    # Load standard spectra if necessary, AB is just calculated analytically so no standard spectrum is required----
    for key, val in filter_dict['standards'].items():
        path = val['path']
        if '.fits' in path:  # If fits file
            with fits.open(path) as hdu:
                standard_df = pd.DataFrame.from_records(hdu[1].data)
            standard_lam, standard_f = standard_df.WAVELENGTH.values, standard_df.FLUX.values
        else:
            standard_txt = np.loadtxt(path)
            standard_lam, standard_f = standard_txt[:, 0], standard_txt[:, 1]
        filter_dict['standards'][key]['lam'] = standard_lam
        filter_dict['standards'][key]['f_lam'] = standard_f

    def ab_standard_flam(l):  # Can just use analytic function for AB spectrum
        f = (const.c.to('AA/s').value / 1e23) * (l ** -2) * 10 ** (-48.6 / 2.5) * 1e23
        return f

    # Load filters------------------------------
    band_weights, zps, offsets = [], [], []
    self.band_dict, self.zp_dict, self.band_lim_dict = {}, {}, {}

    # Prepare NULL band. This is a fake band with a very wide wavelength range used only for padded data points to
    # ensure that these padded data points never fall out of the wavelength coverage of the model. These padded
    # data points do not contribute to the likelihood in any way, this is entirely for computational reasons
    self.band_dict['NULL_BAND'] = 0
    self.zp_dict['NULL_BAND'] = 10  # Arbitrary number
    self.band_lim_dict['NULL_BAND'] = band_wave[0], band_wave[-1]
    band_weights.append(np.ones_like(band_wave))
    zps.append(10)
    offsets.append(0)

    band_ind = 1
    for key, val in filter_dict['filters'].items():
        band, magsys, offset = key, val['magsys'], val['magzero']
        try:
            R = np.loadtxt(val['path'])
        except:
            raise FileNotFoundError(f'Filter response file {val["path"]} not found for {key}')

        # Convert wavelength units if required, model is defined in Angstroms
        units = val.get('lam_unit', 'AA')
        if units.lower() == 'nm':  # Convert from nanometres to Angstroms
            R[:, 0] = R[:, 0] * 10
        elif units.lower() == 'micron':  # Convert from microns to Angstroms
            R[:, 0] = R[:, 0] * 1e4

        band_low_lim = R[np.where(R[:, 1] > 0.01 * R[:, 1].max())[0][0], 0]
        band_up_lim = R[np.where(R[:, 1] > 0.01 * R[:, 1].max())[0][-1], 0]

        # Convolve the bands to match the sampling of the spectrum.
        band_conv_transmission = jnp.interp(band_wave, R[:, 0], R[:, 1], left=0, right=0)
        # band_conv_transmission = scipy.interpolate.interp1d(R[:, 0], R[:, 1], kind='cubic',
        #                                                     fill_value=0, bounds_error=False)(band_wave)

        dlamba = jnp.diff(band_wave)
        dlamba = jnp.r_[dlamba, dlamba[-1]]

        num = band_wave * band_conv_transmission * dlamba
        denom = jnp.sum(num)
        band_weight = num / denom

        band_weights.append(band_weight)

        # Get zero points
        lam = R[:, 0]
        if magsys == 'ab':
            zp = ab_standard_flam(lam)
        else:
            standard = filter_dict['standards'][magsys]
            zp = interp1d(standard['lam'], standard['f_lam'], kind='cubic')(lam)

        int1 = simpson(lam * zp * R[:, 1], x=lam)
        int2 = simpson(lam * R[:, 1], x=lam)
        zp = 2.5 * np.log10(int1 / int2)
        self.band_dict[band] = band_ind
        self.band_lim_dict[band] = [band_low_lim, band_up_lim]
        self.zp_dict[band] = zp
        zps.append(zp)
        offsets.append(offset)
        band_ind += 1

    self.used_band_inds = np.array(list(self.band_dict.values()))
    self.zps = jnp.array(zps)
    self.offsets = jnp.array(offsets)
    self.inv_band_dict = {val: key for key, val in self.band_dict.items()}

    # Get the locations that should be sampled at redshift 0. We can scale these to
    # get the locations at any redshift.
    band_interpolate_locations = jnp.arange(
        0,
        self.spectrum_bins * self.band_oversampling,
        self.band_oversampling
    )

    # Save the variables that we need to do interpolation.
    self.band_interpolate_locations = device_put(band_interpolate_locations)
    self.band_interpolate_spacing = band_spacing
    self.band_interpolate_weights = jnp.array(band_weights)
    self.model_wave = 10 ** model_log_wave
    self.used_band_dict = {val: val for val in self.band_dict.values()}

    self.uv_ind1 = self.model_wave < 2700  # Need to use separate UV term for F99 law below 2700AA
    self.uv_ind2 = (self.model_wave < 2700) & ((1e4 / self.model_wave) >= 5.9)
    self.uv_ind3 = ((1e4 / self.model_wave[self.uv_ind1]) >= 5.9)
    self.uv_x = 1e4 / self.model_wave[self.uv_ind1]

    KD_l = invKD_irr(self.l_knots)
    self.J_l_T = device_put(spline_coeffs_irr(self.model_wave, self.l_knots, KD_l))
    self.KD_t = device_put(invKD_irr(self.tau_knots))
    self._load_hsiao_template()
    self.sim = False  # Keep track of whether data is simulated

    self.ZPT = 27.5  # Zero point
    self.J_l_T = device_put(self.J_l_T)
    self.hsiao_flux = device_put(self.hsiao_flux)
    self.J_l_T_hsiao = device_put(self.J_l_T_hsiao)
    self.xk = jnp.array(
        [0.0, 1e4 / 26500., 1e4 / 12200., 1e4 / 6000., 1e4 / 5470., 1e4 / 4670., 1e4 / 4110., 1e4 / 2700.,
         1e4 / 2600.])
    KD_x = invKD_irr(self.xk)
    self.M_fitz_block = device_put(spline_coeffs_irr(1e4 / self.model_wave, self.xk, KD_x))
