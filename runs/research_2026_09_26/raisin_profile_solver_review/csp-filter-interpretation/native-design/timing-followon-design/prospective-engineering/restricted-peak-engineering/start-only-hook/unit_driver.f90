program unit_start
  implicit none
  integer :: cid=1,it,ip,ierr,check_read,idx,count_records=0
  double precision :: ini(4),snapshot(4),step(4),bounds(2,4),start_value,stored,error,lo,hi,expected,shift
  character(len=20) :: names(4),stored_name
  character(len=64) :: arg,scenario
  call get_command_argument(1,arg)
  read(arg,*) shift
  call get_command_argument(2,scenario)
  names=[character(len=20)::'ITER','ISN','PKMJD','DLMAG']
  do it=1,12
    ini=[dble(it),1d0,57707.80078125d0,42d0]
    step=[0d0,0d0,2d0,.5d0]
    bounds(:,1)=0d0;bounds(:,2)=0d0
    bounds(:,3)=[57671.342999365806d0,57715.067000181196d0]
    bounds(:,4)=[10d0,70d0]
    if(trim(scenario)=='fixed') step(3)=0d0
    if(trim(scenario)=='outside') ini(3)=bounds(2,3)-1d0
    snapshot=ini
    call mninit(5,6,7)
    do ip=1,4
      start_value=ini(ip)
      check_read=0
      if(names(ip)=='PKMJD') then
        call prosp_minuit_peak_start(cid,ini(1),names(1)//char(0),ini(ip),step(ip), &
          bounds(1,ip),bounds(2,ip),start_value,check_read)
      endif
      call mnparm(ip,names(ip),start_value,step(ip),bounds(1,ip),bounds(2,ip),ierr)
      if(ierr/=0)stop 21
      if(check_read==1) then
        call mnpout(ip,stored_name,stored,error,lo,hi,idx)
        ! Deliberate test-only readback corruption; never in the source patch.
        if(trim(scenario)=='missed_readback' .and. it==1)stored=stored-2d0
        if(trim(scenario)=='later_changed' .and. it==2)stored=stored+1d-6
        call prosp_minuit_peak_stored(cid,ini(1),ini(ip),start_value,stored,bounds(1,ip), &
          bounds(2,ip),lo,hi,idx,trim(stored_name)//char(0))
        count_records=count_records+1
      endif
      ! Test-only capture for every parameter, even when hook is disabled.
      call mnpout(ip,stored_name,stored,error,lo,hi,idx)
      expected=ini(ip)
      if(ip==3 .and. it==1)expected=expected+shift
      if(stored/=expected)stop 22
      if(any(ini/=snapshot))stop 23
      if(names(ip)/=stored_name)stop 24
      write(6,'(A,2I3,3ES26.17)') 'UNIT_CAPTURE',it,ip,ini(ip),start_value,stored
    enddo
    ! The same common prior center remains, independent of optimizer start.
    if(ini(3)/=57707.80078125d0)stop 25
  enddo
  if(shift==0d0.and.count_records/=0)stop 26
  if(shift/=0d0.and.count_records/=12)stop 27
  print *, 'UNIT_ALL_PARAMS_PRIOR_CENTER_AND_LATER_ITERATIONS_PASS'
end program
