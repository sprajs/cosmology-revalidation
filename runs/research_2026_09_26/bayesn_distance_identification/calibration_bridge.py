"""Source-linked six-SN photometric-coordinate and passband gate; no fitted outcomes."""
from pathlib import Path
import csv, hashlib, json, re
import numpy as np
from astropy.io import fits
from scipy.integrate import simpson
from scipy.interpolate import interp1d
from ruamel.yaml import YAML

OUT=Path(__file__).resolve().parent
ROOT=OUT.parents[2]
NATIVE=OUT/'official-code/bayesn/bayesn-filters'
RELEASE=ROOT/'sources/repos/djones1040__RAISIN_DataRelease@a383c4b'
yaml=YAML(typ='safe')
MAP={'CSP': {'u':'u_CSP','g':'g_CSP','r':'r_CSP','i':'i_CSP','B':'B_CSP',
             'V':'V_CSP','n':'V_CSP','o':'V_CSP_3009','m':'V_CSP_3014',
             'Y':'Y_RC','y':'Y_WIRC','J':'J_RC1','j':'J_RC2','H':'H_RC','h':'H_WIRC'},
     'PS1MD':{'g':'g_PS1','r':'r_PS1','i':'i_PS1','z':'z_PS1','J':'F125W','H':'F160W'},
     'DES':{'g':'g_DES','r':'r_DES','i':'i_DES','z':'z_DES','J':'F125W','H':'F160W'}}
KCOR={'CSP':'kcor_CSPDR3_BD17.fits','PS1MD':'kcor_PS1MD_NIR.fits','DES':'kcor_DES_NIR.fits'}


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def raw_phot(path):
    head={}; obs=[]
    for line in Path(path).read_text().splitlines():
        a=line.split()
        if not a:continue
        if a[0]=='OBS:':obs.append({'MJD':float(a[1]),'band':a[2],'flux':float(a[4]),'error':float(a[5])})
        elif ':' in a[0]:head[a[0][:-1]]=a[1:]
    return head,obs


def run():
    cfg=yaml.load((NATIVE/'filters.yaml').read_text())
    pilots=list(csv.DictReader((OUT/'pilot-membership.csv').open()))
    used={s:set() for s in MAP}; metadata=[]
    for p in pilots:
        head,obs=raw_phot(ROOT/p['photometry_path'])
        for o in obs:used[p['survey']].add(o['band'])
        metadata.append(dict(CID=p['CID'],N_raw=len(obs),bands=sorted(set(o['band'] for o in obs)),
                             PEAKMJD=float(head['PEAKMJD'][0]),MWEBV=float(head['MWEBV'][0]),
                             zHEL=float(p['zHEL']),zHD=float(p['zHD']),
                             negative_flux_rows=sum(o['flux']<0 for o in obs),
                             nonpositive_error_rows=sum(o['error']<=0 for o in obs),
                             duplicate_peak_headers=sum(x.startswith('PEAKMJD:') for x in (ROOT/p['photometry_path']).read_text().splitlines())))
    (OUT/'release-filters').mkdir(exist_ok=True)
    custom={'standards':{},'filters':{}}
    rows=[]
    for survey,bands in used.items():
        p=RELEASE/'kcor'/KCOR[survey]
        with fits.open(p) as h:
            rt=h['FilterTrans'].data; std=h['PrimarySED'].data; zs=h['ZPoff'].data
            l=np.array(rt.field(0),float)
            for col in std.names[1:]:
                name=f'{survey}_{col}'; path=OUT/'release-filters'/f'{name}_primary.dat'
                np.savetxt(path,np.column_stack([std.field(0),std[col]]),fmt='%.12g')
                custom['standards'][name]={'path':str(path.resolve())}
            for b in sorted(bands):
                oldname=f'CSP-{b}' if survey=='CSP' else (f'{"PS1" if survey=="PS1MD" else "DES"}-{b}/{b}' if b not in 'JH' else f'WFC3_IR_F125W-{b}')
                newname=MAP[survey][b]; new=cfg['filters'][newname]
                zr=zs[np.array([s.strip()==oldname for s in zs['Filter Name']])][0]
                oldT=np.array(rt[oldname],float)
                newR=np.loadtxt(NATIVE/new['path']); ll=newR[:,0]; tt=newR[:,1]
                newT=np.interp(l,ll,tt,left=0,right=0)
                oldnorm=simpson(l*oldT,x=l); newnorm=simpson(l*newT,x=l)
                oldw=l*oldT/oldnorm; neww=l*newT/newnorm
                l1=simpson(np.abs(oldw-neww),x=l)/2
                def primary_flux(magsys,grid):
                    if magsys=='ab':return 2.99792458e18/grid**2*10**(-.4*48.6)
                    info=cfg['standards'][magsys]; sp=NATIVE/info['path']
                    if sp.suffix=='.fits':
                        with fits.open(sp) as sf:s=sf[1].data; return interp1d(s['WAVELENGTH'],s['FLUX'],kind='cubic')(grid)
                    s=np.loadtxt(sp);return interp1d(s[:,0],s[:,1],kind='cubic')(grid)
                oldstd=str(zr['Primary Name']).strip(); oldmag=float(zr['Primary Mag'])
                of=np.interp(l,np.array(std.field(0),float),np.array(std[oldstd],float))
                oldzero=simpson(l*oldT*of,x=l)/oldnorm*10**(.4*oldmag)
                native_f=primary_flux(new['magsys'],ll)
                newzero=simpson(ll*tt*native_f,x=ll)/simpson(ll*tt,x=ll)*10**(.4*new['magzero'])
                scale=oldzero/newzero
                bandkey=f'RAISIN_{survey}_{b}';path=OUT/'release-filters'/f'{bandkey}.dat'
                np.savetxt(path,np.column_stack([l,oldT]),fmt='%.12g')
                custom['filters'][bandkey]={'magsys':f'{survey}_{oldstd}','magzero':oldmag,'path':str(path.resolve())}
                # Native standard integrated through RELEASE curve isolates reference-system change.
                native_on_old=simpson(l*oldT*primary_flux(new['magsys'],l),x=l)/oldnorm*10**(.4*new['magzero'])
                info=dict(survey=survey,raw_band=b,release_filter=oldname,native_filter=newname,custom_filter=bandkey,
                          release_standard=oldstd,release_primary_mag=oldmag,native_standard=new['magsys'],native_primary_mag=new['magzero'],
                          zero_mean_release=oldzero,zero_mean_native=newzero,
                          flux_scale_native_over_release_reference_only=scale,
                          delta_mag_native_minus_release_reference_only=-2.5*np.log10(scale),
                          same_curve_delta_mag=2.5*np.log10(native_on_old/oldzero),
                          photon_weight_half_L1=l1,max_abs_curve_difference=float(np.max(abs(oldT-newT))),
                          release_min_transmission=float(oldT.min()),native_min_transmission=float(tt.min()),
                          release_snphot_zpoff=float(zr['ZPoff(SNpot)']),
                          release_kcor_sha256=sha(p),native_filter_sha256=sha(NATIVE/new['path']),
                          release_curve_sha256=sha(path),native_filter_path=str(NATIVE/new['path']),
                          release_reference_path=str(OUT/'release-filters'/f'{survey}_{oldstd}_primary.dat'),
                          scientific_excluded=(b=='u'))
                rows.append(info)
    with (OUT/'calibration-bridge.csv').open('w') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    outyaml=YAML();outyaml.dump(custom,(OUT/'release-filter-config.yaml').open('w'))
    (OUT/'pilot-photometry-metadata.json').write_text(json.dumps(metadata,indent=2))
    (OUT/'calibration-bridge-result.json').write_text(json.dumps(dict(
        status='No photometry transformed or posterior fitted; native and released reference systems audited separately.',
        filter_count=len(rows),pilot_metadata=metadata,
        nonzero_SNphot_offsets=[r for r in rows if r['release_snphot_zpoff']!=0],
        reference_delta_range_mag=[min(r['delta_mag_native_minus_release_reference_only'] for r in rows),max(r['delta_mag_native_minus_release_reference_only'] for r in rows)],
        shape_difference_range_half_L1=[min(r['photon_weight_half_L1'] for r in rows),max(r['photon_weight_half_L1'] for r in rows)],
        flux_error_rule='If coordinate scaling is used, error_new=S*error_old and C_new=S*C_old*S.T; differing filters additionally require SED transport.',
        KCOR_primary_units='per Angstrom: checked current pinned SNANA kcor.c wr_fits_PRIMARY plus source spectral magnitude values; ignore obsolete per10A header description',
        forward_bridge='Use exact released KCOR curves, primary spectra and primary magnitudes as custom BayeSN filters; this preserves raw FLUXCAL/error coordinates and avoids any unsupported scalar transport.',
        caveat='Controlled release-native bridge, not execution-linked 2024 BayeSN preprocessing. Rounded10Angstrom KCOR tables require convergence check.'),indent=2))
    print(json.dumps({r['survey']+'_'+r['raw_band']:{'dm_mmag':1000*r['delta_mag_native_minus_release_reference_only'], 'halfL1':r['photon_weight_half_L1']} for r in rows},indent=2))


if __name__=='__main__':run()
