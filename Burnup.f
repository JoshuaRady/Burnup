!---------------------------------------------------------------------------------------------------
! Burnup.f
! Burnup Source Code, Original Fixed Form
!
! Original code by: Frank A. Albini
! Edited and annotated by: Joshua M. Rady
! Woodwell Climate Research Center
! 9/2023
!
! The following is annotated source code for the original version of the Burnup model from:
!
! Program BURNUP, a simulation model of the burning of large woody natural fuels.
! Frank A. Albini
! Missoula, MT: USDA Forest Service. Unpublished report. Research Grant INT-92754-GR. 1994.
! Appendix B.
!
! The original report is held by the National Forest Service Library in Fort Collins.
!
! Code Transcription:
! 	A scanned version of the original report in PDF format was used as the basis for this code as
! the original source code could not be found in digital form.  The optical character recognition
! (OCR) for the PDF appendix holding the code has very extensive errors so I have edited it
! line-by-line and tested compilation and output.  Some missing pages we added to the scan after
! correspondence with the librarians.
!
! 	The original Burnup source code is self-described as being in Fortran IV, though there are
! numerous features from Fortran 66 and 77.  It is in fixed form form.  I have used this to help
! interpret the whitespace in the original text.  I have tried to preserve even arbitrary whitespace
! for fidelity with the print version but in some places it is approximate.  Some lines are longer
! than the 72 character limit imposed for fixed form.  For these line I commented out the original
! code line and provided a modified line that follows below.
!
! 	I have used tabs for code indenting, where indenting is clearing indicated in the original code,
! and for alignment of comments.  This is to aid with modernization and porting of this code.  While
! tabs are not valid Fortran characters, compilers tend to tolerate them.  I have used 4-space
! equivalent tabs for layout purposes and to reduce line lengths but the original code uses ~3-4
! space indenting.  Lines that use tabs to maintain lengths under 72 characters but look longer due
! to tab alignment have been marked with JMR_Length_OK_Tabs.
!
! 	The original code uses the Fortran IV comment character 'c' at the first line position.  In most
! places that is cut off in the scan and I figured it out by context.  I have replaced the 'c'
! comment character with '!c' for original comments, again to aid in updating and porting.  Comments 
! starting with just '!' have been added by me.  Page breaks and page numbers from the PDF have been
! marked in comments to make comparison to the orignal code easier.
!
! Caveats:
! 	This code compiles, runs, and produces output but transcription errors may remain, some of which
! may affect results.
!
! 	There is no licence provided for the original code.  It is though to be open by provenance, but
! that man not be correct.
!---------------------------------------------------------------------------------------------------


! Pg. 75: Start of source code.


      program BURNUP
      parameter ( MAXNO = 10 )
      parameter ( MAXKL = MAXNO * ( MAXNO + 1 ) / 2 + MAXNO )
      parameter ( MXSTEP = 20 )

      real*4 wdry( maxno ) , ash( maxno ) , htval( maxno )
      real*4 fmois( maxno ) , dendry( maxno ) , sigma( maxno )
      real*4 cheat( maxno ) , condry( maxno ) , alfa( maxno )
      real*4 tpig( maxno ) , tchar( maxno )
      real*4 flit( maxno ) , fout( maxno ) , work( maxno )
      real*4 elam( maxno , maxno ) , alone( maxno )
      real*4 area( maxno ) , fint( maxno )
      real*4 xmat( maxkl ) , tdry( maxkl )
      real*4 tign( maxkl ) , tout( maxkl )
      real*4 wo( maxkl ) , wodot( maxkl )
      real*4 diam( maxkl ) , ddot( maxkl )
      real*4 qcum( maxkl ) , tcum( maxkl )
      real*4 acum( maxkl ) , qdot( maxkl , mxstep )
      integer key( maxno)
      character*12 parts( maxno ) , list( maxno ) , infile , outfil
      logical nohist

      fimin = 0.1
      nruns = 0

10    call GETDAT( infile , outfil , parts , wdry , ash , htval ,
     +             fmois , dendry , sigma , cheat , condry , tpig ,
     +             tchar , number , maxno , fi , ti , u , d , tpamb ,
     +             ak, r0 , dr , dt , ntimes , wdf , dfm, nruns ,
     +             area, fint )

11    write ( * , 12 )
12    format(' Enter 1 to store fire history, 0 to skip ' ,$)
      read( * , * ) ihist
      nohist = ( ihist .EQ. 0 )
      if ( ( ihist .NE. 1 ) .AND. ( .NOT. nohist ) ) goto 11

      call ARRAYS( maxno, number , wdry, ash, dendry , fmois ,
     +             sigma, htval , cheat , condry, tpig , tchar ,
     +             diam, key , work , ak , elam , alone , xmat , wo ,
     +             maxkl , parts , list , area )

      now = 1
      tis = ti
      call DUFBRN( wdf , dfm , dfi , tdf )

      call START( tis , mxstep , now , maxno , number , wo , alfa ,
     +            dendry , fmois , cheat , condry , diam , tpig ,
     +            tchar , xmat , tpamb , tpdry , fi , flit , fout ,
     +            tdry, tign , tout , qcum , tcum , acum , qdot ,
     +            ddot, wodot , work , u , d , r0 , dr , ch2o ,
     +            ncalls , maxkl )

      if( tis .LT. tdf ) then


! -- Pagebreak --
! Pg. 76:


      fid = dfi
      else
      fid = 0.

      end if
      if( nohist ) goto 14

      call STASH( tis , now, maxno, number , outfil , fi , flit ,
     + fout, wo , wodot , diam , ddot , tdry , tign,
     + tout , fmois , maxkl , nun )

14    call FIRINT( wodot , ash , htval , maxno , number , maxkl , area ,
     +             fint, fi )

      if( fi .LE. fimin ) goto 22

20    call STEP( dt , mxstep , now , maxno , number , wo , alfa ,
     +           dendry , fmois , cheat , condry , diam , tpig ,
     +           tchar , xmat , tpamb , tpdry , fi , flit , fout ,
     +           tdry , tign , tout , qcum , tcum , acum , qdot ,
     +           ddot , wodot , work , u , d , r0 , dr , ch2o ,
     +           ncalls , maxkl , tis , fint , fid )

      now = now + 1
      tis = tis + dt
      if ( tis .LT. tdf ) then
      	fid = dfi
      else
      	fid = 0.
      end if

      call FIRINT( wodot , ash , htval , maxno , number , maxkl , area ,
     + fint, fi )

      if ( fi .LE. fimin ) goto 22

      if ( nohist ) goto 21

      call STASH( tis, now, maxno, number, outfil, fi, flit,
     +            fout , wo , wodot , diam , ddot , tdry , tign,
     +            tout , fmois , maxkl , nun)

21    continue

      if ( now .LT. ntimes ) goto 20

22    close( nun )

      outfil = 'SUMMARY.DAT'
      call SUMMARY( outfil , number , maxno , maxkl , parts , nun,
     +tis , ak , wdry , fmois , sigma , tign , tout , xmat , wo , diam )
      close( nun )
25    write(*,30)
30    format(' Exercise completed. Do another = 1 , quit = 0 ',$)
      read (*,*,err=25) next


! -- Pagebreak --
! Pg. 77:


      if( next .EQ. 0) then
      	close( 66)
      	stop' Terminated'
      end if
      if ( next . NE. 1 ) goto 25
      close( 66 )
      goto 10

      end


      subroutine DUFBRN( wdf , dfm , dfi , tdf )

!c Duff burning rate (ergo, intensity) and duration

      dfi = 0.
      tdf = 0.
      if( ( wdf .LE. 0. ) .OR. ( dfm .GE. 1.96 ) ) return
      dfi = 11.25 - 4.05 * dfm
      ff = 0.837 - 0.426 * dfm
      tdf = 1.e+04 * ff * wdf / ( 7.5 - 2.7 * dfm )
      return
      end

! -- Pagebreak --
! Pg. 78:


      subroutine GETDAT( infile , outfil , parts , wdry , ash , htval ,
!     +                   fmois , dendry , sigma , cheat , condry , tpig , ! JMR_NOTE: Too long!
!     +                   tchar , number , maxno , fi , ti , u , d , tpamb , ! JMR_NOTE: Too long!
!     +                   ak , r0 , dr , dt , ntimes , wdf , dfm , nruns , ! JMR_NOTE: Too long!
     +                  fmois , dendry , sigma , cheat , condry , tpig , ! JMR_MOD
     +                tchar , number , maxno , fi , ti , u , d , tpamb , ! JMR_MOD
     +                  ak , r0 , dr , dt , ntimes , wdf , dfm , nruns , ! JMR_MOD
     +                   area , fint)

!c This routine prompts user for input data in the form of file names
!c or the direct input of numerical quantities.
!c Routine makes consistency checks, tests nominal quantities, and
!c allows review / revision of data before starting program,
!c and allows archiving of run conditions for future use.

      real*4 wdry( maxno ) , ash( maxno ) , htval( maxno )
      real*4 fmois( maxno ) , dendry( maxno ) , sigma( maxno )
      real*4 cheat( maxno ) , condry( maxno )
      real*4 tpig( maxno ) , tchar( maxno )
      real*4 area( maxno ) , fint( maxno )
      character*12 parts( maxno ) , infile , outfil , name1
      logical test

      data small / 1.e-06 / , big / 1.e+06 /
      data ash1 / 0.0001 / , ash2 / 0.1 /
      data htv1 / 1.0e+07 / , htv2 / 3.0e+07 /
      data fms1 / 0.01 / , fms2 / 3.0 /
      data den1 / 200.0 / , den2 / 1000.0 /
      data sig1 / 4.0 / , sig2 / 1.0e+04 /
      data cht1 / 1000. 0 / , cht2 / 2000. 0 /
      data con1 / 0.025 / , con2 / 0.25 /
      data tig1 / 200.0 / , tig2 / 400.0 /
      data tch1 / 250.0 / , tch2 / 500.0 /
      data fir1 / 100. / , fir2 / 1.0e+05 /
      data ti1 / 10.0 / , ti2 / 200.0 /
      data u1 / 0.0 / , u2 / 5.0 /
      data d1 / 0.1 / , d2 / 2.0 /
      data tam1 / -40.0 / , tam2 / 40.0 /
      data wdf1 / 0.1 / , wdf2 / 30.0 /
      data dfm1 / 0.1 / , dfm2 / 1.972 /


      if( nruns .GT. 0 ) goto 01
      do i = 1 , maxno
      	wdry( i ) = 0.
      	ash( i ) = 0.
      	htval( i ) = 0.
      	fmois( i ) = 0.
      	dendry( i ) = 0.
      	sigma( i ) = 0.
      	cheat( i ) = 0.
      	condry( i ) = 0.
      	tpig( i ) = 0.
      	tchar( i ) = 0.
      	area( i ) = 0.
      	fint( i ) = 0.
      	parts( i ) = 'No data for'
      end do
      number = 0


! -- Pagebreak --
! Pg. 79:


      fi = 0.
      ti = 0.
      u = 0.
      d = 0.
      tpamb = 0.
      ak = 0.
      r0 = 0.
      dr = 0.
      dt = 0.
      ntimes = 0
      wdf = 0.
      dfm = 2.

01    write(*,02)
02    format(' Enter 1 for keyboard entry of data'/
     +       ' Enter 2 to fetch data stored in files'/
     +       ' Enter 0 to terminate program now ' , $)
      read(*,*,err=01) in
      if( in .EQ. 0) stop' Program ended'
      if( in .EQ. 2 ) goto 1000
      if( in .NE. 1 ) goto 01

      nruns = nruns + 1

03    write(*,04) maxno
04    format(' How many [ no more than'i3' ] fuel components ?    ',$)
      read(*,*,err=03) numb
      if( ( numb .LE. 0 ) .OR. ( numb .GT. maxno ) ) goto 03
      number = numb

      nn = 0

0505  nn = nn + 1
05    write(*,06) nn
06    format(//' Fuel component number'i3 ' properties list'/
     +         ' Name [ 12 characters or fewer ]    ', $ )
      read(*,08) name1
      parts( nn ) = name1
08    format(a12)

09    write(*,10) name1
10    format(1x,a12 ' dry loading, kg / sq m ',$)
      read(*,*,err=09) drywt
      call retry( drywt , small , big , test )
      if( test) goto 09
      wdry( nn) = drywt

11    write(*,12) name1
12    format(1x,a12 ' total mineral ash fraction    ',$)
      read(*,*,err=11) ashes
      call retry( ashes , ash1 , ash2 , test)
      if( test ) goto 11
      ash( nn) = ashes

13    write(*,14) name1


! -- Pagebreak --
! Pg. 80:


14    format(1x,a12' heat of combustion, J / kg    ',$)
      read(*,*,err=13) hots
      call retry( hots , htv1 , htv2 , test)
      if( test) goto 13
      htval( nn ) = hots

15    write(*,16) name1
16    format(1x,a12' moisture content, fraction    ',$)
      read(*, *,err=15) fms
      call retry( fms , fms1 , fms2 , test )
      if( test) goto 15
      fmois( nn ) = fms


17    write(*, 18) name1
18    format(1x,a12' ovendry mass density, kg / cubic m    ',$)
      read(*,*,err=17) dryd
      call retry( dryd , den1 , den2 , test )
      if( test) goto 17
      dendry( nn ) = dryd

19    write(*,20) name1
20    format(1x,a12' surface / volume ratio , inverse m    ',$)
      read(*,*,err=19) sigs
      call retry( sigs , sig1 , sig2 , test )
      if( test ) goto 19
      sigma( nn ) = sigs

21    write(*,22) name1
! 22    format(1x,a12' specific heat capacity ovendry , J / kg deg K    ',$) ! JMR_NOTE: Too long!
22    format(1x,a12' specific heat capacity ovendry , J / kg deg K  ',$) ! JMR_MOD
      read(*,*,err=21) cpd
      call retry( cpd , cht1 , cht2, test )
      if( test ) goto 21
      cheat( nn ) = cpd

23    write(*,24) name1
24    format(1x,a12' ovendry thermal conductivity , W / m deg K    ',$)
      read(*,*,err=23) cond
      call retry( cond , con1 , con2 , test )
      if( test) goto 23
      condry( nn ) = cond

25    write(*,26) name1
26    format(1x,a12' ignition temperature, deg C    ',$)
      read(*,*,err=25) tigi
      call retry( tigi , tig1 , tig2 , test )
      if( test) goto 25
      tpig( nn ) = tigi + 273.0

27    write(*,28) name1
28    format(1x,a12' char [end pyrolysis] temperature, deg C    ',$)
      read (*,*,err=27) tchi
      call retry( tchi , tch1 , tch2 , test )
      if( test ) goto 27
      tchar( nn ) = tchi + 273.0


! -- Pagebreak --
! Pg. 81:


      if( nn .LT. number ) goto 0505
      if( ( next .EQ. 1 ) .AND. ( ido .EQ. 1 ) ) goto 50
      
29    write(*,30)
30    format(' Enter igniting surface fire intensity , kW / sq m   ',$)
      read(*,*,err=29) fir
      call retry( fir , fir1 , fir2 , test )
      if( test ) goto 29
      fi = fir

31    write (*,32)
32    format(' Enter igniting surface fire residence time , s   ',$)
      read(*,*,err=31) tir
      call retry( tir , ti1 , ti2 , test )
      if( test ) goto 31
      ti = tir

33    write(*,34)
34    format(' Windspeed at top of fuelbed , m / s   ',$)
      read(*,*,err=33) uin
      call retry( uin , u1 , u2 , test )
      if ( test ) goto 33
      u = uin

35    write(*,36)
36    format(' Depth of fuelbed , m   ',$)
      read(*,*,err=35) din
      call retry( din , d1 , d2 , test )
      if( test ) goto 35
      d = din

37    write(*,38)
38    format(' Ambient temperature , deg C   ',$ )
      read(*,*,err=37) tamb
      call retry( tamb , tam1 , tam2 , test )
      if( test ) goto 37
      tpamb = tamb + 273.0

39    write(*,40)
40    format(' Dimensionless area influence factor [ak parameter]   ',$)
      read(*,*,err=39) ak

41    write(*,42)
42    format(' Fire environment minimum temperature parameter r0   ',$)
      read(*,*,err=41) r0

43    write(*,44)
!44    format(' Fire environment increment temperature parameter dr   ',$) ! JMR_NOTE: Too long!
44    format(' Fire environment increment temperature parameter dr  ',$) ! JMR_MOD
      read(*,*,err=43) dr

45    write(*,46)
46    format(' Time step for integration of burning rates , s   ',$)
      read(*,*,err=45) dt

47    write(*,48)


! -- Pagebreak --
! Pg. 82:


48    format(' Number time steps [ must exceed 1 ]   ',$)
      read(*,*,err=47) ntimes
      if( ntimes .LE. 1 ) goto 47

49    write(*,491)
491   format(' Duff dry weight loading, kg / sq m   ',$)
      read(*,*,err=49) wdf
      call retry( wdf , wdf1 , wdf2 , test )
      if( test ) goto 49

492   write(*, 493)
493   format(' Duff moisture fraction   ',$)
      read(*,*,err= 492) dfm
      call retry( dfm , dfm1 , dfm2 , test )
      if( test) goto 492

!c All data collected --- now time to review and revise as needed

50    write(*,52)
52    format(' All data can be reviewed and revised now; Enter:'/
     +       ' 0 to archive and/or execute present data'/
     +       ' 1 to review fuel component properties'/
     +       ' 2 to review igniting fire and environmental data'/
     +       ' 3 to review internal and control variables'/
     +       ' 4 to terminate program now   ',$)
      read(*,*,err=50) next
      if( next .EQ. 4) stop' Terminated'
      if( next .EQ. 0) goto 2000
      if( ( next .LT. 1 ) .OR. ( next .GT. 3 ) ) goto 50
      if( next .EQ. 1 ) then

!c Review fuel component properties

53    	write(*,54)
54    	format(' Enter 0 to terminate review & revision of fuel data'/
     +         ' Enter 1 to select a fuel component for review'/
     +         ' Enter 2 to add or delete a fuel component   ',$)
      	read(*,*,err=53) ido
      	if( ido .EQ. 0 ) goto 50
      	if( ido .EQ. 2 ) goto 66
      	if( ido .NE. 1 ) goto 53
55    	write(*,56) ( n , parts( n ) , wdry( n ) , n = 1 , number )
56    	format(i3,3x,a12,'   loading = 'e12.3'   kg / sq m')
57    	write(*,58)
58    	format(' Number component to be revised [ 0 = prev menu ]   ',$)
      	read(*,*,err=57) nn
      	if( nn .EQ. 0 )goto 53
      	if( ( nn .LT. 1 ) .OR. ( nn .GT. number ) ) goto 55
59    	write(*,60) parts( nn ) , wdry( nn ) , ash( nn ) ,
     +	        htval( nn ) , fmois( nn ) , dendry( nn ) ,
     +	        sigma( nn ) , cheat( nn ) , condry( nn ) ,
     +	        ( tpig( nn ) - 273.0 ) , ( tchar( nn ) - 273.0 )
60    	format(' index parameter data entered'/
     +         '   1'10x'name'4x,a12 /
     +         '   2'7x'loading  ',e12.3/


! -- Pagebreak --
! Pg. 83:


     +         '    3'7x'ash frac ',e12.3/
     +         '    4'7x'heat comb',e12.3/
     +         '    5'7x'mois frac',e12.3/
     +         '    6'7x'density  ',e12.3/
     +         '    7'7x'surf/vol ',e12.3/
     +         '    8'7x'heat capy',e12.3/
     +         '    9'7x'conduct y',e12.3/
     +         '   10'7x'ig temp C',e12.3/
     +         '   11'7x'char temp',e12.3/
     +         '   enter number param to change [ 0 = prev menu ]   ',$)
      	read(*,*,err=59) idn
      	if( idn .EQ. 0 ) goto 55
      	if( idn .GT. 11 ) goto 59
      	if( idn .EQ. 1 ) then
      		write(*,*)' Enter new name - 12 char or fewer'
      		read(*,08) parts( nn )
      		goto 59
      	end if
61    	write(*,6 2)
62    	format(' New value = ?   '$)
      	read(*,*,err=61) value
      	if( idn .EQ. 2 ) wdry( nn ) = value
      	if( idn .EQ. 3 ) ash( nn) = value
      	if( idn .EQ. 4 ) htval( nn ) = value
      	if( idn .EQ. 5 ) fmois( nn ) = value
      	if( idn .EQ. 6 ) dendry( nn) = value
      	if( idn .EQ. 7 ) sigma( nn) = value
      	if( idn .EQ. 8 ) cheat( nn) = value
      	if( idn .EQ. 9 ) condry( nn ) = value
      	if( idn .EQ. 10 ) tpig( nn ) = value + 273.0
      	if( idn .EQ. 11 ) tchar( nn	) = value + 273.0
      	goto 59

!c Add or delete a fuel component

66    	write(*,68)
68    	format(' Enter - 1 to delete a fuel component'/
     +         ' Enter 0 if no more additions or deletions'/
     +         ' Enter + 1 to add a fuel component   ', $)
      	read(*,*,err=66) ido
      	if ( ido .EQ. 0 ) goto 50
      	if( ido .EQ. -1 ) then

!c Delete a fuel component

69    		write(*,56) ( n , parts( n ) , wdry( n ) , n = 1 , number )
      		write(*,70)
70    		format(' Index number of component to delete   ', $)
      		read(*,*,err=69) ind
      		if( ( ind .LE. 0 ) .OR. ( ind .GT. number ) ) goto 66
      		if( ind .EQ. number ) then
      			number = number - 1
      			goto 66
      		end if
      		parts( ind ) = parts( number )


! -- Pagebreak --
! Pg. 84:


      		wdry( ind ) = wdry( number )
      		ash( ind ) = ash( number )
      		htval( ind ) = htval( number )
      		fmois( ind ) = fmois( number )
      		dendry( ind ) = dendry( number )
      		sigma( ind ) = sigma( number )
      		cheat( ind ) = cheat( number )
      		condry( ind ) = condry( number )
      		tpig( ind ) = tpig( number )
      		tchar( ind ) = tchar( number )
      		number = number - 1
      		goto 66
      	end if

!c Add a fue1 component

      	if( ido .NE. 1 ) goto 66
      	nn = number
      	number = number + 1
      	goto 0505
      end if

      if ( next .EQ. 2 ) then

!c Review igniting fire and environmental data

71    	write(*,72) fi
72    	format('  1 - Igniting fire intensity, kW/ sq m      'e12.3)
      	write(*,73) ti
73    	format('  2 - Igniting fire residence time, s        'e12.3)
      	write (*,74) u
74    	format('  3 - windspeed at top of fuel bed, m / s    'e12.3)
      	write(*,75) d
75    	format('  4 - fuelbed depth, m                       'e12.3)
      	write(*,76) ( tpamb - 273.0 )
76    	format('  5 - ambient temperature, deg C             'e12.3)
77    	write(*,78)
78    	format(' Index of parameter to change [ 0 = prev menu ]   ',$)
      	read(*,*,err=77) ind
      	if( ind .EQ. 0 ) goto 50
      	if( ( ind .LT. 1 ) .OR. ( ind .GT. 5 ) ) goto 71
      	write(*,62)
      	read(*,*,err=71) value
      	if( ind .EQ. 1 ) fi = value
      	if( ind .EQ. 2 ) ti = value
      	if( ind .EQ. 3 ) u = value
      	if( ind .EQ. 4 ) d = value
      	if( ind .EQ. 5 ) tpamb = value + 273.0
      	goto 71
      end if

!c input quantity next .EQ. 3  by fall through to this point, so
!c Review internal and conttol variables

79    write(*,80) ak


! -- Pagebreak --
! Pg. 85:


      write(*,81) r0
      write(*,82) dr
      write(*,83) dt
      write(*,84) wdf
      write(*,85) dfm
      write(*,86) ntimes
      write(*,87)
80    format('  1 - Area influence factor [ ak parameter ]       'e12.3)
81    format('  2 - Fire environment temperature parameter r0    'e12.3)
82    format('  3 - Fire environment temperature parameter dr    'e12.3)
83    format('  4 - Time step to integrate burning rates , s     'e12.3)
84    format('  5 - Duff dry weight loading , kg / sq m          'e12.3)
85    format('  6 - Duff moisture content , fraction dry weight  'e12.3)
86    format('  7 - Number time steps                            'i4)
87    format('  Index of parameter to change [ 0 = prev menu ]   ',$)
      read(*,*,err=79 ) ind
      if( ind .EQ. 0 ) goto 50
      if( ( ind .LT. 1 ) .OR. ( ind .GT. 7 ) ) goto 79
      write(*,62)
      if( ind .LT. 7) read(*,*,err=79) value
      if( ind .EQ. 1 ) ak = value
      if( ind .EQ. 2 ) r0 = value
      if( ind .EQ. 3 ) dr = value
      if( ind .EQ. 4 ) dt = value
      if( ind .EQ. 5 ) wdf = value
      if( ind .EQ. 6 ) dfm = value
      if( ind .EQ. 7) read(*,*,err=79) ntimes
      goto 79

!c End of keyboard entry and review-data / revise-data segment
!c Following are file retrieval [ 1000 ] and archiving [ 2000 ] segments

1000  write(*, 1001)
1001  format(' Enter name of file holding fuel component data   ',$)
      read(*,08,err=1000) infile
1002  open(77,file=infile,status='OLD',form='FORMATTED',err=1022)
      read(77,1003) number
1003  format(i6)
      do n = 1 , number
      	read(77,2004) parts( n )
      	read(77,2005) wdry( n ) , ash( n ) , htval( n ) , fmois( n ) ,
     +	               dendry( n ) , sigma( n ) , cheat( n ) ,
     +	               condry( n ) , tpig( n ) , tchar( n )
      end do
      close(77)
      write(*,1004)
1004  format(' Enter name of file holding igniting fire, environmental'/
     +       ', and program control parameters   ',$)
      read(*,08) infile
1005  open(77,file=infile,status='OLD' ,form='FORMATTED' ,err=1025)
      read(77,2008) fi , ti , u , d , tpamb
      read(77,2009) ak , r0 , dr , dt , ntimes
      read(77,2010,end=1006) wdf , dfm
      close(77)
      goto 50


! -- Pagebreak --
! Pg. 86:


1006  wdf = 0.
      dfm = 0.
      close(77)
      goto 50
1022  write(*,1030) infile
1030  format(' Error opening input file  ', a12'  reenter file name')
      read(*,08,err=1022) infile
      goto 1002
1025  write(*,1030) infile
      read(*,08,err=1025) infile
      goto 1005

!c Archive current run conditions?

2000  write(*,2100)
2100  format(' Enter 2 for previous menu'/
     +       ' Enter 1 to archive current data set'/
     +       ' Enter 0 to execute without archiving   '$)
      read(*,*,err=2000) ido
      if( ido .EQ. 2 ) goto 50
      if( ido .EQ. 0 ) return
      if( ido .NE. 1 ) goto 2000
2001  write(*,2002)
2002  format(' Enter file name [ 12 char or fewer ] for storage of the'/
     +       ' fuel component data currently in use   ',$)
      read(*,08,err=2001) outfil
2022  open(66,file=outfil,status='UNKNOWN',form='FORMATTED',err= 2021)
      write(66,2003) number
2003  format(i6,' Fuel component properties -- 3 records per entry')
      do n= 1 , number
      	write(66,2004) parts( n )
      	write(66, 2005) wdry( n ) , ash( n ) , htval( n ) , fmois( n ) ,
     +	                dendry( n ) , sigma( n ) , cheat( n ) ,
     +	                condry( n ) , tpig( n ) , tchar( n )
2004	format(2x,a12)
2005	format(2x,5e15.3/2x,5e15.3)
      end do
      close(66)
2006  write(*,2007)
2007  format(' Enter file name [ 12 char or fewer ] for storage of the'/
!     +       ' igniting fire, environmental, and control data used   ',$ ! JMR_NOTE: Too long!
     +       ' igniting fire, environmental, and control data used  ',$) ! JMR_MOD
      read(*,08,err=2006) outfil
2027  open(66,file=outfil,status='UNKNOWN',form='FORMATTED',err=2026)
      write(66,2008) fi , ti , u, d, tpamb
2008  format(2x,5e15.3)
      write(66,2009) ak , r0 , dr , dt , ntimes
2009  format(2x,4e15.3, i6)
      write(66,2010) wdf , dfm
2010  format(2x,2e15.3)
      close(66)
      goto 2000
2021  write(*,2030) outfil
2030  format(' Error opening output file  'a12'  enter another name')
      goto 2001
2026  write(*,2030) outfil


! -- Pagebreak --
! Pg. 87:


      read(*,08,err=2026) outfil
      goto 2027

      end



      subroutine RETRY( value , valo , vahi , test )
      logical test

      test = ( value .GT. valo ) .AND. ( value .LT. vahi )
      if( test ) then
      	test = .FALSE.
      	return
      end if

      test = ( value .LE. valo )
      if( test ) then
01    	write(*,*)' Value entered seems small'
      	write(*,*)' Enter 0 to proceed with small value'
      	write(*,*)' Enter 1 to provide a new value'
      	read(*,*,err=01 ) in
      	if( in .EQ. 0) then
      		test = .FALSE.
      		return
      	end if
      	if( in .NE. 1 ) goto 01
      	return
      end if

      test = ( value .GE. vahi )
      if( test ) then
02    	write(*,*)' Value entered seems large'
      	write(*,*)' Enter 0 to proceed with large value'
      	write(*,*)' Enter 1 to provide a new value'
      	read(*,*,err=02) in
      	if( in .EQ. 0) then
      		test = .FALSE.
      		return
      	end if
      	if( in .NE. 1) goto 02
      	return
      end if
      end


! -- Pagebreak --
! Pg. 88:  ARRAYS() starts here.


      subroutine ARRAYS( maxno , number , wdry , ash , dendry , fmois ,
     +                   sigma , htval , cheat , condry , tpig , tchar ,
     +                   diam , key , work , ak , elam , alone , xmat ,
     +                   wo , maxkl , parts , list , area )

!c Orders the fuel description arrays according to the paradigm described in
!c subroutine SORTER and computes the interaction matrix xmat from the array
!c elam and the list alone returned from subroutine OVLAPS.  Parameters in
!c arrays are defined in terms of initial values as:
!c		wdry		ovendry mass loading , kg / sq m
!c		ash			mineral content , fraction dry mass
!c		dendry		ovendry mass density , kg / cu m
!c		fmois		moisture content , fraction dry mass
!c		sigma		surface to volume ratio , 1 / m
!c		htval		low heat of combustion , J / kg
!c		cheat		specific heat capacity , ( J / K) / kg dry mass
!c		condry		thermal conductivity , W / m K , ovendry
!c		tpig		ignition temperature , K
!c		tchar		char temperature , K
!c		diam		initial diameter , m [ by interaction pairs ]
!c		key			ordered index list
!c		work		workspace array
!c		elam		interaction matrix from OVLAPS
!c		alone		noninteraction fraction list from OVLAPS
!c		xmat		consolidated interaction matrix
!c		wo			initial dry loading by interaction pairs
!c		area		fraction of site area expected to be covered at
!c					least once by initial planform area of ea size


      character*12 parts( maxno ) , list( maxno )
      real*4 wdry( maxno ) , ash( maxno ) , dendry( maxno )
      real*4 fmois( maxno ) , sigma( maxno ) , htval( maxno )
      real*4 cheat( maxno ) , condry( maxno ) , tpig( maxno)
      real*4 tchar( maxno ) , work( maxno )
      real*4 elam( maxno , maxno) , alone ( maxno )
      real*4 area( maxno )
      real*4 xmat( maxkl ) , diam( maxkl ) , wo( maxkl )

      integer key( maxno )

      loc( k , l ) = k * ( k + 1 ) / 2 + l

      call SORTER( maxno , number , sigma, fmois , dendry , key)

      do j = 1 , number
      	k = key( j )
      	list( j ) = parts( k )
      end do
      do j = 1 , number
      	parts( j ) = list ( j )
      end do
      do j = 1 , number


! -- Pagebreak --
! Pg. 89:


      	k = key( j )
      	work( j ) = wdry( k )
      end do
      do j = 1 , number
      	wdry( j ) = work( j )
      end do

      do j = 1 , number
      	k = key( j )
      	work( j ) = ash( k )
      end do
      do j = 1 , number
      	ash( j ) = work( j )
      end do

      do j = 1 , number
      	k = key( j )
      	work( j ) = htval( k )
      end do
      do j = 1 , number
      	htval( j ) = work( j )
      end do

      do j = 1 , number
      	k = key( j )
      	work( j ) = cheat( k )
      end do
      do j = 1 , number
      	cheat( j ) = work( j )
      end do

      do j = 1 , number
      	k = key( j )
      	work( j ) = condry( k )
      end do
      do j = 1 , number
      	condry( j ) = work( j )
      end do

      do j = 1 , number
      	k = key( j )
      	work( j ) = tpig( k )
      end do
      do j = 1 , number
      	tpig( j ) =  work (j)
      end do

      do j = 1 , number
      	k = key( j )
      	work( j ) = tchar( k )
      end do
      do j = 1 , number
      	tchar( j ) = work( j )
      end do

! -- Pagebreak --
! Pg. 90:


      call OVLAPS( wdry , sigma , dendry , ak , number , maxno , maxkl ,
     +                                   xmat , elam , alone , area )

      do k = 1 , number
      	diak = 4. / sigma( k )
      	wtk = wdry( k )
      	kl = loc( k , 0 )
      	diam( kl ) = diak
      	xmat( kl ) = alone( k )
      	wo( kl ) = wtk * xmat( kl )
      	do j = 1 , k
      		kj = loc( k , j )
      		diam( kj ) = diak
      		xmat( kj ) = elam( k , j )
      		wo( kj ) = wtk * xmat( kj )
      	end do
      end do

      return
      end


! -- Pagebreak --
! Pg. 91:


      subroutine SORTER( maxno , number , sigma , fmois , dryden , key )

!c Sorts fuel element list in order of increasing size (decreasing sigma)
!c For elements with same size order determined on increasing moisture
!c content (fmois). If items have same size and moisture content, order
!c on the basis of increasing mass density (dryden). "number" elements are
!c included in the list, which has a maximum length of "maxno". The integer
!c list: key( j ) , j = 1 , number holds the indices in order, so other
!c fuel parameters can be ordered and associated as necessary.

      real*4 sigma( maxno ) , fmois( maxno ) , dryden( maxno )
      integer key( maxno )
      logical diam , mois , dens , tied

      do j = 1 , maxno
      	key ( j ) = j
      end do

!c Replacement sort: order on increasing size , moisture , density

      do j = 2, number
      	s = 1. / sigma( j )
      	fm = fmois( j )
      	de = dryden( j )
      	keep = key( j )
      	do i = ( j - 1 ) , 1 , -1
      		usi = 1. / sigma( i )
      		diam = ( usi .LT. s )
      		if( diam ) goto 10
      		tied = ( usi .EQ. s )
      		if( .NOT. tied ) goto 05
      		mois = ( fmois( i ) .LT. fm )
      		if( mois ) goto 10
      		tied = ( fmois( i ) .EQ. fm )
      		if( .NOT. tied ) goto 05
      		dens = ( dryden ( i ) .LE. de )
      		if( dens ) goto 10
05    		sigma( i + 1 ) = sigma( i )
      		fmois( i + 1) = fmois( i )
      		dryden( i + 1 ) = dryden( i )
      		key( i + 1 ) = key( i )
      	end do
      	i = 0
10    	sigma( i + 1 ) = 1. / s
      	fmois( i + 1 ) = fm
      	dryden( i + 1 ) = de
      	key( i + 1 ) = keep
      end do

      return
      end

! -- Pagebreak --
! Pg. 92:  This page did not have OCR applied.  I extracted the page and ran OCR on it myself.


      subroutine OVLAPS( dryld , sigma , dryden, ak , number , maxno ,
     +                              maxkl , beta , elam , alone , area )
!c Computes the interaction matrix elam( j , k) which apportions the
!c influence of smaller and equal size pieces on eaqh size class for the
!c purpose of establishing the rates at which the elemnts burn out.
!c Input quantities are dryld, the ovendry mass per unit area of each
!c element available for burning after the passage of the igniting surface
!c fire; sigma, the particle's surface/ volume ratio , and dryden, the
!c ovendry mass density of the particle; ak a dimensionless parameter that
!c scales the planform area of a particle to its area of influence. There
!c are "number" separate particle classes, of a maximum number = maxno.
!c It is assumed that the 1ist are ordered on size class (nonincreasing
!c surface/ volume ratio). List "alone" gives the fraction of each loading
!c that is not influenced by any other category.

      real*4 dryld( maxno ) , sigma( maxno ) , dryden( maxno)
      real*4 beta( maxkl ) , elam( maxno , maxno) , alone( maxno )
      real*4 area( maxno )

      loc( k , l ) = k * ( k + 1 ) / 2 + l

      pi = abs( acos( -1. ) )
      do j = 1 , number
      	alone( j ) = 0.
      	do k = 1 , j
      		kj = loc( j , k)
      		beta( kj ) = 0.
      	end do
      	do k = 1 , number
      		elam( j , k ) = 0.
      	end do
      end do

      do k = 1 , number
      	siga = ak * sigma( k ) / pi
      	do l = 1 , k
      		kl = loc( k , l )
      		a = siga * dryld( l ) / dryden( l )
      		if ( k .EQ. 1 ) then
      			bb = 1. - exp( - a )
      			area( k ) = bb
      		end if
      		if ( k .NE. 1 ) bb = min( 1. , a )
      		beta( kl ) = bb
      	end do
      end do

      if ( number .EQ. 1 ) then
      	elam( 1 , 1 ) = beta( 2 )
      	alone( 1 ) = 1. - elam( 1 , 1 )
      	return
      end if

      do k = 1 , number


! -- Pagebreak --
! Pg. 93:


      	frac = 0.
      	do l = 1 , k
      		kl = loc( k , l )
      		frac = frac + beta( kl )
      	end do
      	if ( frac .GT. 1. ) then
      		do l = 1 , k
      			kl = loc( k , l )
      			elam( k , l ) = beta( kl ) / frac
      		end do
      		alone( k ) = 0.
      	else
      		do l = 1 , k
      			kl = loc( k, l )
      			elam( k , l ) = beta( kl )
      		end do
      		alone( k ) = 1. - frac
      	end if
      end do

      return
      end


! -- Pagebreak --
! Pg. 94:


      subroutine START( dt , mxstep , now , maxno , number , wo , alfa ,
     +                 dendry , fmois , cheat , condry , diam , tpig ,
     +                 tchar , xmat , tpamb , tpdry , fi , flit , fout ,
     +                 tdry , tign , tout , qcum , tcum , acum , qdot ,
     +                 ddot , wodot , work , u , d , r0 , dr , ch2o ,
     +                 ncalls , maxkl )

!c This routine initializes variables prior to starting sequence of calls
!c to subroutine STEP.  On input here, fi is area intensity of spreading
!c fire , dt is the residence time for the spreading fire.  Only physical
!c parameters specified are the fuel array descriptors. To call STEP ,
!c one must initialize the following variables.

!c Input parameters:
!c	dt =		spreading fire residence time , sec
!c	mxstep =	max dimension of historical sequences
!c	now = 		index marks end of time step
!c	maxno =		max number of fuel components
!c	number =	actual number of fuel components
!c	wo =		current ovendry loading for the larger of
!c				each component pair, kg / sq m
!c	alfa =		dry thermal diffusivity of component , sq m / s
!c	dendry =	ovendry density of component , kg /cu m
!c	fmois =		moisture fraction of component
!c	cheat =		specific heat capacity of component , J / kg K
!c	condry =	ovendry thermal conductivity , W / sq m K
!c	diam =		current diameter of the larger of each
!c				fuel component pair , m
!c	tpig =		ignition temperature ( K) , by component
!c	tchar =		end - pyrolysis temperature ( K ) , by component
!c	xmat =		table-of-influence fractions between components
!c	tpamb =		ambient temperature ( K)
!c	fi =		current fire intensity ( site avg )  , kW / sq m

!c Parameters updated [input and output]
!c	ncalls = counter of calls to this routine
!c                   = 0 on first call or reset
!c                   cumulates after first call
!c	flit =		fraction of each component currently alight
!c	fout =		fraction of each component currently gone out
!c	tdry =		time of drying start of the larger of each
!c				fuel component pair
!c	tign =		ignition time for the larger of each
!c				fuel component pair
!c	tout =		burnout time of larger component of pairs
!c	qcum =		cumulative heat input to larger of pair , J / sq m
!c	tcum =		cumulative temp integral for qcum ( drying )
!c	acum =		heat pulse area for historical rate averaging
!c	qdot =		history ( post ignite ) of heat transfer rate
!c				to the larger of each component pair
!c	ddot =		diameter reduction rate , larger of pair , m / s
!c	wodot =		dry loading loss rate for larger of pair


! -- Pagebreak --
! Pg. 95:


!c Constant parameters
!c	u = 		mean horizontal windspeed at top of fuelbed
!c	d =			fuelbed depth
!c	r0 = 		minimum value of mixing parameter
!c	dr =		max - min value of mixing parameter
!c	ch2o =		specific heat capacity of water, J / kg K
!c	hvap =		heat of vaporization of water J / kg
!c	tpdry =		temperature ( all components ) start drying ( K )

      real*4 wo( maxkl ) , alfa( maxno ) , dendry( maxno )
      real*4 fmois( maxno ) , cheat( maxno )
      real*4 condry( maxno ) , diam( maxkl ), tpig( maxno)
      real*4 tchar( maxno ) , xmat( maxkl )
      real*4 tdry( maxkl ) , tign( maxkl )
      real*4 tout( maxkl ) , qcum( maxkl )
      real*4 tcum( maxkl ) , acum( maxkl )
      real*4 qdot( maxkl , mxstep ) , ddot( maxkl)
      real*4 wodot( maxkl ) , flit( maxno ) , fout( maxno )
      real*4 work( maxno )

      loc( k , l ) = k * ( k + 1 ) / 2 + l

      data rindef / 1.e+30 /

      ch2o = 4186.
      tpdry = 353.

!c Initialize time varying quantities and, set up work( k )
!c The diameter reduction rate of fuel component k is given
!c by the product of the rate of heat tranefer to it, per
!c unit surface area, and the quantity work( k )

      do k = 1 , number
      	fout( k ) = 0.
      	flit( k ) = 0.
      	alfa( k ) = condry( k) / ( dendry( k) * cheat( k ) )
!c		effect of moisture content on burning rate (scale factor)
      	delm = 1.67 * fmois( k )
!c		effect of component mass density (empirical)
      	heatk = dendry( k ) / 446.
!c		empirical burn rate factor, J / cu m - K
      	heatk = heatk * 2.01e+06 * ( 1. + delm )
!c		normalize out driving temperature difference (Tfire - Tchar)
!c		to average value of lab experiments used to find above constants
      	work( k ) = 1. / ( 255. * heatk )
      	do l = 0 , k
      		kl = loc( k , l )
      		tout( kl ) = rindef
      		tign( kl ) = rindef
      		tdry( kl ) = rindef
      		tcum( kl ) = 0.
      		qcum( kl ) = 0.
      	end do
      end do


! -- Pagebreak --
! Pg. 96: This page did not have OCR applied.  I extracted the page and ran OCR on it myself.


!c Make first estimate of drying start times for all components
!c These times are usually brief and make little or no difference

      r = r0 + 0.25 * dr
      tf = tempf( fi , r , tpamb )
      ts = tpamb
      if( tf .LE. ( tpdry + 10. ) ) stop' Igniting fire cannot dry fuel'
      thd = ( tpdry - ts ) / ( tf - ts )
      tx = 0.5 * (ts + tpdry)

      do k = 1 , number
      	factor = dendry( k ) * fmois( k )
      	conwet = condry( k ) + 4.27e-04 * factor
      	do l = 0 , k
      		kl = loc( k , l )
      		dia = diam( kl )
      		call heatx( u , d , dia , tf , tx , hf , hb , conwet , en )
      		call DRYTIM( en , thd , dryt )
      		cpwet = cheat( k ) + fmois( k ) * ch2o
      		fac = ( ( 0.5 * dia ) ** 2 ) / conwet
      		fac = fac * dendry( k ) * cpwet
      		dryt = fac * dryt
      		tdry( kl ) = dryt
      	end do
      end do

!c Next, determine which components are alight in spreading fire

      tsd = tpdry

      do k = 1 , number
      	c = condry( k )
      	tigk = tpig( k )
      	do l = 0 , k
      		kl = loc( k , l )
      		dryt = tdry( kl )
      		if( dryt .GE. dt) goto 10
      		dia = diam( kl )
      		ts = 0.5*( tsd + tigk )
      		call heatx( u , d , dia , tf , ts , hf , hb , c , e )
      		tcum( kl ) = max( ( tf - ts ) * ( dt - dryt ) , 0. )
      		qcum( kl ) = hb * tcum( kl )
      		if ( tf .LE. ( tigk + 10. ) ) goto 10
      		call TIGNIT( tpamb , tpdry , tpig ( k ) , tf , condry ( k ) , ! JMR_Length_OK_Tabs
     +		             cheat( k ) , fmois( k ) , dendry( k ) , hb ,
     +		             dtign )
      		trt = dryt + dtign
      		tign( kl ) = 0.5 * trt
      		if( dt .GT. trt ) flit( k ) = flit( k ) + xmat( kl )
10    		continue
      	end do
      end do

      nlit = 0


! -- Pagebreak --
! Pg. 97:


      trt = rindef

!c Determine minimum ignition time and verify ignition exists

      do k = 1 , number
      	if( flit( k ) .GT. 0. ) nlit = nlit + 1
      	do l = 0 , k
      		kl = loc( k , l )
      		trt = min( trt , tign( kl ) )
      	end do
      end do

      	if( nlit .EQ. 0 ) stop' START ignites no fuel'

!c Deduct trt from all time estimates , resetting time origin

      do k = 1 , number
      	do l = 0 , k
      		kl = loc( k , l )
      		if( tdry( kl ) .LT. rindef ) then
      			tdry( kl ) = tdry ( kl ) - trt
      		end if
      		if( tign( kl ) .LT. rindef) then
      			tign( kl ) = tign( kl ) - trt
      		end if
      	end do
      end do

!c Now go through all component pairs and establish burning rates
!c for all the components that are ignited; extrapolate to end time dt

      do k = 1 , number
      	if( flit( k ) .EQ. 0. ) then
      		do l = 0 , k
      			kl = loc( k , l )
      			ddot( kl ) = 0.
      			tout( kl ) = rindef
      			wodot( kl ) = 0.
      		end do
      	else
      		ts = tchar( k )
      		c = condry( k )
      		do l = 0 , k
      			kl = loc( k , l )
      			dia = diam( kl )
      			call heatx( u , d , dia , tf , ts , hf , hb , c , e )
      			qdot( kl , now ) = hb * max( ( tf - ts ) , 0. )
      			aint = ( c / hb) ** 2
      			ddt = dt - tign( kl )
      			acum( kl ) = aint * ddt
      			ddot( kl ) = qdot( kl , now ) * work( k )
      			tout( kl ) = dia / ddot( kl )
      			dnext = max( 0. , ( dia - ddt * ddot( kl ) ) )
      			wnext = wo( kl ) * ( ( dnext / dia ) ** 2 )
      			wodot( kl ) = ( wo( kl ) - wnext ) / ddt


! -- Pagebreak --
! Pg. 98:


      			diam( kl ) = dnext
      			wo( kl ) = wnext
      			df = 0.
      			if ( dnext .LE. 0. ) then
      				df = xmat( kl )
      				wodot( kl ) = 0.
      				ddot( kl ) = 0.
      			end if
      			flit( k ) = flit( k ) - df
      			fout( k ) = fout( k ) + df
      		end do
      	end if
      end do


      ncalls = 0

      return
      end


! -- Pagebreak --
! Pg. 99:


      subroutine FIRINT( wodot , ash , htval , maxno , number , maxkl ,
     +                   area , fint , fi )
!c Computes fi = site avg fire intensity given the burning rates of all
!c interacting pairs of fuel components [ wodot ] , the mineral ash content
!c of each component [ ash ] , the heat of combustion value [ htval ] for
!c each , and the number of fuel components [ number ] , where max - maxno.
!c fi is in kW / sq m, while htval is in J / kg.

!c fint( k ) is the correction to fi to adjust
!c the intensity level to be the local value where size k is burning.

      real*4 ash( maxno ) , htval( maxno )
      real*4 wodot( maxkl ) , area( maxno ) , fint( maxno )
      data small / 1.e-06 /

      loc( k , l ) = k * ( k + 1 ) / 2 + l

      sum = 0.
      do k = 1 , number
      	wdotk = 0.
      	do l = 0, k
      		kl = loc( k , l )
      		wdotk = wdotk + wodot( kl )
      	end do
      	term = ( 1. - ash( k ) ) * htval( k ) * wdotk * 1.e-03
      	ark = area( k )
      	if ( ark .GT. small ) then
      		fint( k ) = term / ark - term
      	else
      		fint( k ) = 0.
      	end if
      	sum = sum + term
      end do

      fi = sum
      return
      end

! -- Pagebreak --
! Pg. 100:


      subroutine STASH( time , now , maxno , number , outfil , fi ,
     +                  flit , fout , wo , wodot , diam , ddot ,
     +                  tdry , tign , tout , fmois , maxkl , nun)

!c This routine stashes output from the BURNUP model package on a snapshot
!c basis.  Every time it is called, it "dumps" a picture of the status of
!c each component of the fuel complex, as a table of interacting pairs.

      real*4 wo( maxkl ) , wodot ( maxkl )
      real*4 diam( maxkl ) , ddot( maxkl )
      real*4 flit( maxno ) , fout( maxno)
      real*4 tdry( maxkl) , tign( maxkl )
      real*4 tout( maxkl ) , fmois( maxno)
      character*12 outfil, histry
      logical snaps

      loc( k , l ) = k * ( k + 1 ) / 2 + l

      if( now .EQ. 1) then
      	wd = 0.
      	wg = 0.
      	do m = 1 , number
      		fmm = fmois( m )
      		wdm = 0.
      		do n = 0 , m
      			mn = loc( m, n)
      			wdm = wdm + wo( mn )
      		end do
      		wgm = wdm * ( 1. + fmm )
      		wd = wd + wdm
      		wg = wg + wgm
      	end do
      	wd0 = wd
      	wg0 = wg
      	snaps = .FALSE.
      	nun = 77
      	mum = 66 
      	histry = 'HISTORY.DAT'
      	open( mum, file=histry , status='UNKNOWN' , form='FORMATTED' )
1000  	write( * , 1001 )
1001  	format(' Enter 1 to stash fuelbed snapshots, O to skip ',$)
      	read( * , * ) istash
      	if( istash .EQ. 0 ) goto 06
      	if( istash .NE. 1 ) goto 1000
      	snaps = .TRUE.
01    	write( * , 02 )
02    	format(' File name [ max 12 char ] for STASH output    ', $)
      	read(*,03) outfil
03    	format (a12)
      	open(nun, file=outfil, status='NEW', form='FORMATTED' , err=04)
      	goto 06
04    	write(*,05) outfil
05    	format(' Error opening file    'a12/
     +         ' Enter 0 to try a new name'/
     +         ' Enter 1 to overwrite existing file'/


! -- Pagebreak --
! Pg. 101:


     +         ' Enter 2 to terminate program now    ',$)
      	read(*,*) ido
      	if( ido .EQ. 2 ) stop' Program ended'
      	if( ido .EQ. 0 ) goto 01
      	if( ido .NE. 1 ) goto 04
      	open(nun,file=outfil,status='UNKNOWN',form='FORMATTED',err=04)
      end if

06    continue

      if( snaps ) write( nun , 07 ) time , fi
07    format( 8e10.3)

      wd = 0.
      wg = 0.
      do m = 1 , number
      	fmm = fmois( m )
      	wdm = 0.
      	if( snaps ) write( nun , 07 ) time , flit( m ) , fout( m )
      	do n = 0 , m
      		mn = loc( m , n )
      		wdm = wdm + wo( mn )
      		if( snaps ) write( nun, 07 ) time , wo( mn ) ,
     +                  wodot( mn ) , diam( mn ) , ddot( mn ) ,
     +                   tdry( mn ) , tign( mn ) , tout( mn )
      	end do
      	wgm = wdm * ( 1. + fmm )
      	wd = wd + wdm
      	wg = wg + wgm
      end do
      wgf = wg / wg0
      wdf = wd / wd0

      write( mum , 07 ) time , wg , wd , wgf , wdf , fi

      return
      end


! -- Pagebreak --
! Pg. 102: This page did not have OCR applied.  I extracted the page and ran OCR on it myself.


      subroutine TIGNIT( tpam, tpdr , tpig, tpfi , cond ,
     +                   chtd, fmof , dend, hbar , tmig )

!c tpam = ambient temperature , K
!c tpdr = fuel temperature at start of drying , K
!c tpig = fuel surface temperature at iinition, K
!c tpfi = fire enviroriment temperature , K
!c cond = fuel ovendry thermal conductivity, W / m K
!c chtd = fuel ovendry specific heat capacity, J / kg K
!c fmof = fuel moisture content , fraction dry weight
!c dend = fuel ovendry density, kg / cu m
!c hbar = effective film heat transfer coefficient [< HEATX] W / sq m K
!c tmig = predicted time to piloted ignition , s

      data a03 / -1.3371565 / , a13 / 0.4653628 / , a23 / -0.1282064 /
      data pinv / 2.125534 / , small / 1.e-06 /
      data hvap / 2.177e+06 / , cpm / 4186. / , conc / 4.27e-04 /

!c approximate function of beta to be solved is ff( x ) where
!c x  =  1 / ( 1 + p * beta )  { Hastings, Approximations for
!c digital computers } and we employ      pinv  =  1 / p

      ff( x) = b03 + x * ( a13 + x * ( a23 + x) )

!c radiant heating equivalent form gives end condition fixes beta

      b03 = a03 * ( tpfi - tpig ) / ( tpfi - tpam )

!c find x that solves ff( x) = 0 ; method is binary search

      xlo = 0.
      xhi = 1.
01    xav = 0.5 * ( xlo + xhi )
      fav = ff( xav )
      if( abs( fav ) .LE. small ) goto 10
      if( fav .LT. 0. ) xlo = xav
      if( fav .GT. 0. ) xhi = xav
      goto 01

10    beta = pinv * ( 1. - xav ) / xav
      conw = cond + conc * dend * fmof
      dtb = tpdr - tpam
      dti = tpig - tpam
      ratio = ( hvap + cpm * dtb ) / ( chtd * dti )
      rhoc = dend * chtd * ( 1. + fmof * ratio )
      tmig = ( ( beta / hbar ) ** 2 ) * conw * rhoc
      
      return
      end


! -- Pagebreak --
! Pg. 103:


      subroutine DRYTIM( enu , theta , tau )

!c Given a Nusselt number ( enu , actually Biot number = h D / k)
!c and the dimensionless temperature rise required for the start
!c of surface drying ( theta ), returns the dimerisionless time ( tau )
!c needed to achieve it. The time is given multiplied by thermal
!c diffusivity and divided by radius squared. Solution by binary search.

      data a / 0.7478556 / , b / 0.4653628 / , c / 0.1282064 /
      data p / 0.47047 /
      
      f( h ) = h * ( b - h * ( c - h ) ) - ( 1. - theta ) / a

      xl = 0.
      xh = 1.

      do n = 1 , 15
      	xm = 0.5 * ( xl + xh )
      	if( f( xm ) .LT. 0. ) then
      		xl = xm
      	else
      		xh = xm
      	end if
      end do

      x = ( 1. / xm - 1. ) / p
      tau = ( 0.5 * x / enu ) **2
      return
      end
      

! -- Pagebreak --
! Pg. 104:


      subroutine HEATX( u , d , dia , tf , ts , hfm , hbar , cond , en )

!c Given horizontal windspeed u at height d [top of fuelbed], cylindrical
!c fuel particle diameter dia, fire environment temperature tf, and mean
!c surface temperature, ts, subroutine returns film heat transfer coefficient
!c hfm and an "effective" film heat transfer coefficient including radiation
!c heat transfer, hbar.  Using the wood's thermal conductivity, cond, the
!c modified Nusselt number [ en ] used to estimate onset of surface drying
!c is returned as well.

      data g / 9.8 /
      data vis / 7.5 e - 05 /
      data a / 8.75 e - 03 /
      data b / 5.75 e - 05 /
      data rad / 5.67 e - 08 /
      data fmfac / 0.382 /
      data hradf / 0.5 /

      hfm = 0.
      if( dia .LE. b ) goto 10
      v = sqrt( u * u + 0.53 * g * d )
      re = v * dia / vis
      enuair = 0.344 * ( re ** 0.56 )
      conair = a + b * tf
      fac = sqrt( abs( tf - ts ) / dia )
      hfmin = fmfac * sqrt( fac )
      hfm = max( ( enuair * conair / dia ) , hfmin )
10    hrad = hradf * rad * ( tf + ts ) * ( tf * tf + ts * ts )
      hbar = hfm + hrad
      en = hbar * dia / cond
      return
      end


! -- Pagebreak --
! Pg. 105:


      function TEMPF( q , r , tamb )
!c Returns a fire environment temperature , TEMPF , given the fire intensity
!c q in kW / square meter , the ambient temperature tamb in Kelvins, and the
!c dimensionless mixing parameter r.

      data err / 1. e - 04 / , aa / 20. /
      term = r / ( aa * q )
      rlast = r
10    den = 1. + term * ( rlast + 1. ) * ( rlast * rlast + 1. )
      rnext = 0.5 * ( rlast + 1. + r / den )
      if( abs( rnext - rlast) .LT. err ) then
      	tempf = rnext * tamb
      	return
      end if
      rlast = rnext
      goto 10
      end


! -- Pagebreak --
! Pg. 106:


      subroutine STEP( dt , mxstep , now , maxno , number , wo , alfa ,
     +                 dendry , fmois , cheat , condry , diam , tpig ,
     +                 tchar , xmat , tpamb , tpdry , fi , flit , fout ,
     +                 tdry , tign , tout , qcum , tcum , acum , qdot ,
     +                 ddot , wodot , work , u , d , r0 , dr , ch2o ,
     +                 ncalls , maxkl , tin , fint , fid )

!c Updates status of all fuel component pairs and returns a snapshot

!c Input parameters:

!c	tin =			start of current time step
!c	dt =			time step , sec
!c	mxstep =		max dimension of historical seouences
!c	now =			index marks end of time step
!c	maxno =			max number of fuel components
!c	number =		actual number of fuel components
!c	wo =			current ovendry loading for the larger of
!c					each component pair, kg / sq m
!c	alfa =			dry thermal diffusivity of component , sq m / s
!c	dendry =		ovendry density of component , kg / cu m
!c	fmois =			moisture fraction of component
!c	cheat =			specific heat capacity of component , J / kg K
!c	condry =		ovendry thermal conductivity , w / sq m K
!c	diam =			current diameter of the larger of each
!c					fuel component pair , m
!c	tpig =			ignition temperature ( K ) , by component
!c	tchar =			end - pyrolysis temperature ( K) , by component
!c	xmat =			table of influence fractions between components
!c	tpamb =			ambient temperature ( K )
!c	tpdry =			temperature ( all components ) start drying ( K )
!c	fi =			current fire intensity ( site avg ) , kW / sq m
!c	work( k ) =		factor of heat transfer rate hbar * (Tfire - Tebar)
!c					that yields ddot ( k )
!c	fint ( k ) =	correction to fi to compute local intensity
!c					that may be different due to k burning
!c	fid =			fire intensity due to duff burning ... this is
!c					used to up the fire intensity for fuel pieces
!c					that are burning without interacting with others
!c	plus the following constants and bookk keeping parameters
!c	u, d, r0 , dr , ch20 , ncalls , maxkl

!c Parameters updated [input and output)

!c	ncalls =		counter of calls to this routine ...
!c							= 0 on first call or reset
!c							cumulates after first call
!c	flit =			fraction of each component currently alight
!c	fout =			fraction of each component currently gone out
!c	tdry =			time of drying start of the larger of each
!c					fuel component pair
!c	tign =			ignition time for the larger of each
!c					fuel component pair
!c	tout =			burnout time of larger component of pairs
!c	qcum =			cumulative heat input to larger of pair , J / sq m


! -- Pagebreak --
! Pg. 107:


!c	tcum =			cumulative temp integral for qcum ( drying )
!c	acum =			heat pulse area for historical rate averaging
!c	qdot =			history ( post ignite ) of heat transfer rate
!c					to the larger of component pair , W / sq m
!c	ddot =			diameter reduction rate , larger of pair , m / s
!c	wodot =			dry loading loss rate for larger of pair

!c Constant parameters
!c	u =				mean horizontal windspeed at top of fuelbed
!c	d =				fuelbed depth
!c	r0 =			minimum value of mixing parameter
!c	dr =			max - min value of mixing parameter
!c	ch2o =			specific heat capacity of water , J / kg K
!c	hvap =			heat of vaporization of water , J / kg

      real*4 wo( maxkl ) , alfa( maxno ) , dendry( maxno )
      real*4 fmois( maxno ) , cheat( maxno )
      real*4 condry( maxno ) , diam( maxkl ) , tpig( maxno )
      real*4 tchar( maxno ) , xmat( maxkl )
      real*4 tdry( maxkl ) , tign( maxkl )
      real*4 tout( maxkl ) , qcum( maxkl )
      real*4 tcum( maxkl ) , acum( maxkl )
      real*4 qdot( maxkl , mxstep ) , ddot( maxkl )
      real*4 wodot( maxkl) , flit( maxno ) , fout( maxno)
      real*4 work( maxno ) , fint( maxno )
      real*4 rindef / 1.e+30 /

      logical flag

      loc( k , l ) = k * ( k + 1 ) / 2 + l

      ncalls = ncalls + 1
      tnow = tin
      tnext = tnow + dt
!c tifi = time when fire ignition phase ended ( at now = 1 )
      tifi = tnow - float( now - 1 ) * dt
      next = now + 1

      do k = 1 , number
      	c = condry( k)
      	do l = 0 , k
      		kl = loc( k , l )
      		tdun = tout( kl )

!c See if k of ( k , l ) pair burned out

      		if ( tnow .GE. tdun ) then
      			ddot( kl) = 0.
      			wodot( kl ) = 0.
      			goto 10
      		end if
      		if( tnext .GE. tdun) then
      			tgo = tdun - tnow
      			ddot( kl ) = diam( kl ) / tgo


! -- Pagebreak --
! Pg. 108:


      			wodot( kl ) = wo( kl ) / tgo
      			wo( kl ) = 0.
      			diam( kl ) = 0.
      			goto 10
      		end if

!c k has not yet burned out ... see if k of ( k, l ) pair is ignited

      		tlit = tign( kl )
      		if( tnow .GE. tlit ) then
      			ts = tchar( k )
      			if( l .EQ. 0 ) then
      				r = r0 + 0.5 * dr
      				gi = fi + fid
      			end if
      			if ( ( l .NE. 0 ) .AND. ( l .NE. k ) ) then
      				r = r0 + 0.5 * ( 1. + flit( l ) ) * dr
      				gi = fi + fint( k ) + flit( l ) * fint( l )
      			end if
      			if( l .EQ. k ) then
      				r = r0 + 0.5 * ( 1. + flit( k ) ) * dr
      				gi = fi + flit ( k ) * fint( k )
      			end if
      			tf = tempf( gi , r , tpamb )
      			dia = diam( kl )
      			call heatx( u, d , dia , tf , ts , hf , hb , c , e )
      			qqq = hb * max( ( tf - ts ) ,  0. )
      			tst = max( tlit , tifi )
      			nspan = max( l , nint( ( tnext - tst ) / dt ) )
      			if( nspan .LE. mxstep) qdot( kl , nspan ) = qqq
      			if( nspan .GT. mxstep) then
      				do mu = 2, mxstep
      					qdot( kl , mu - 1 ) = qdot ( kl , mu )
      				end do
      				qdot( kl , mxstep) = qqq
      			end if
      			aint = ( c / hb ) ** 2
      			acum( kl ) = acum( kl ) + aint * dt
      			tav1 = tnext - tlit
      			tav2 = acum( kl ) / alfa( k )
      			tav3 = ( ( dia / 4. ) ** 2 ) / alfa ( k )
      			tavg = min( tav1 , tav2 , tav3 )
      			index = 1 + min( nspan , mxstep )
      			qdsum = 0.
      			tspan = 0.
      			deltim = dt
01    			index = index - 1
      			if( index .EQ. 1 ) deltim = tnext - tspan - tlit
      			if( ( tspan + deltim) .GE. tavg ) deltim = tavg - tspan
      			qdsum = qdsum + qdot( kl , index ) * deltim
      			tspan = tspan + deltim
      			if( ( tspan .LT. tavg ) .AND. ( index .GT. 1 ) ) goto 01
      			qdavg = max ( qdsum / tspan , 0. )
      			ddot( kl ) = qdavg * work( k )
      			dnext = max( 0. , dia - dt * ddot( kl ) )


! -- Pagebreak --
! Pg. 109:


      			wnext = wo( kl ) * ( ( dnext / dia) ** 2 )
      			if( ( dnext .EQ. 0. ) .AND. ( ddot( kl ) .GT. 0. ) ) then ! JMR_Length_OK_Tabs
      				tout(kl ) = tnow + dia / ddot( kl )
      			end if
      			if( ( dnext .GT. 0. ) .AND. ( dnext .LT. dia ) ) then
      				rate = dia / ( dia - dnext)
      				tout( kl ) = tnow + rate * dt
      			end if
      			if( qdavg .LE. 20. ) tout( kl ) = 0.5*( tnow + tnext )
      			ddt = min( dt , ( tout( kl ) - tnow ) )
      			wodot( kl ) = ( wo( kl ) - wnext ) / ddt
      			diam( kl ) = dnext
      			wo( kl ) = wnext
      			goto 10
      		end if

!c See if k of ( k, l ) has reached outer surface drying stage yet

      		dryt = tdry( kl )
      		if( ( tnow .GE. dryt) .AND. ( tnow .LT. tlit) ) then
      			if( l .EQ. 0 ) then
      				r = r0
      				gi = fi + fid
      			end if
      			if( l .EQ. k ) then
      				r = r0
      				gi = fi
      			end if
      			if( ( l .NE. 0 ) .AND. ( l .NE. k) ) then
      				r = r0 + 0.5 * flit( l ) * dr
      				gi = fi + flit ( l ) * fint( l )
      			end if
      			tf = tempf( gi, r , tpamb )
      			ts = tpamb
      			dia = diam ( kl )
      			call heatx( u , d , dia , tf , ts , hf , hb , c , e )
      			dtemp = max( 0. , ( tf - ts ) )
      			dqdt = hb * dtemp
      			qcum( kl ) = qcum( kl ) + dqdt * dt
      			tcum( kl ) = tcum( kl ) + dtemp * dt
      			dteff = tcum( kl ) / ( tnext - dryt )
      			heff = qcum( kl ) / tcum( kl )
      			tfe = ts + dteff
      			dtlite = rindef
      			if( tfe .LE. ( tpig ( k ) + 10. ) ) goto 02
      			call TIGNIT( tpamb , tpdry , tpig( k ) , tfe ,
     +           condry( k ) , cheat( k ) , fmois( k ) , dendry( k ) ,
     +           heff , dtlite )
02    			tign( kl ) = 0.5 * ( dryt + dtlite )

!c If k will ignite before time step over , must interpolate

      			if( tnext .GT. tign( kl ) ) then
      				ts = tchar( k )
      				call heatx( u , d , dia , tf , ts , hf , hb , c , e ) ! JMR_Length_OK_Tabs

! -- Pagebreak --
! Pg. 110:


      				qdot( kl , 1 ) = hb * max( ( tf - ts ) , 0. )
      				qd = qdot( kl , 1 )
      				ddot( kl ) = qd * work( k )
      				delt = tnext - tign( kl )
      				dnext = max( 0. , dia - delt * ddot( kl ) )
      				wnext = wo( kl ) * ( ( dnext / dia ) ** 2 )
      				if( dnext .EQ. 0. ) tout( kl ) = tnow + dia / ddot(kl) ! JMR_Length_OK_Tabs
      				if( ( dnext .GT. 0. ) .AND. ( dnext .LT. dia ) ) then
      					rate = dia / ( dia - dnext )
      					tout( kl ) = tnow + rate * dt
      				end if
      				if ( tout( kl ) .GT.  now ) then
      					ddt = min( dt , ( tout( kl ) - tnow ) )
      					wodot( kl ) = ( wo( kl ) - wnext ) / ddt
      				else
      					wodot( kl ) = 0.
      				end if
      				diam( kl ) = dnext
      				wo( kl ) = wnext
      			end if
      			goto 10
      		end if

!c If k of ( k , l ) still coming up to drying temperature , accumulate
!c heat input and driving temperature difference , predict drying start

      		if( tnow .LT. dryt) then
      			factor = fmois( k ) * dendry( k )
      			conwet = condry( k ) + 4.27e-04 * factor
      			if( l .EQ. 0 ) then
      				r = r0
      				gi = fi + fid
      			end if
      			if( l .EQ. k ) then
      				r = r0
      				gi = fi
      			end if
      			if( ( l .NE. 0 ) .AND. ( l .NE. k )) then
      				r = r0 + 0.5 * flit( l ) * dr
      				gi = fi + flit( l ) * fint( l )
      			end if
      			tf = tempf( gi , r , tpamb )
      			if( tf .LE. ( tpdry + 10. )) goto 10
      			dia = diam( kl )
      			ts = 0.5*( tpamb + tpdry)
      			call heatx( u, d, dia , tf , ts , hf , hb , c , e )
      			dtcum = max( ( tf - ts ) * dt , 0. )
      			tcum( kl ) = tcum( kl ) + dtcum
      			qcum( kl ) = qcum( kl ) + hb * dtcum
      			he = qcum( kl ) / tcum( kl )
      			dtef = tcum( kl ) / tnext
      			thd = ( tpdry - tpamb ) / dtef
      			if( thd .GT. 0.9 ) goto 10
      			biot = he * dia / conwet
      			call DRYTIM( biot , thd , dryt )


! -- Pagebreak --
! Pg. 111:


      			cpwet = cheat( k ) + ch2o * fmois( k )
      			fac = ( ( 0.5 * dia ) ** 2 ) / conwet
      			fac = fac * cpwet * dendry( k )
      			tdry( kl ) = fac * dryt
      			if( tdry( kl ) .LT. tnext ) then
      				ts = tpdry
      				call heatx( u , d , dia , tf , ts , hf , hb , c , e ) ! JMR_Length_OK_Tabs.
      				dqdt = hb * ( tf - ts )
      				delt = tnext - tdry( kl )
      				qcum( kl ) = dqdt * delt
      				tcum( kl ) = ( tf - ts ) * delt
      				tbar = 0.5 * ( tpdry + tpig( k ) )

!c See if ignition to occur before time step complete

      				if( tf .LE. ( tpig( k ) + 10. ) ) goto 10
      				call TIGNIT( tpamb , tpdry , tpig( k ) , tf ,
     +		         condry( k ) , cheat( k ) , fmois( k ) , dendry( k ) , ! JMR_Length_OK_Tabs
     +               hb , dtlite )
      				tign( kl ) = 0.5 * ( tdry( kl ) + dtlite )
      				if( tnext .GT. tign( kl ) ) then
      					ts = tchar( k)
      					qdot( kl , 1 ) = hb * max ( ( tf - ts ) , 0. )
      				end if
      			end if
      		end if
10    		continue
      	end do
      end do

!c Update fractions ignited and burned out , to apply at next step start

      do k = 1 , number
      	flit ( k ) = 0.
      	fout ( k ) = 0.
      	do l = 0 , k
      		kl = loc( k , l )
      		flag = ( tnext .GE. tign( kl ) )
      		if ( flag .AND. ( tnext .LE. tout( kl ) ) ) then
      			flit( k) = flit( k ) + xmat( kl )
      		end if
      		if ( tnext .GT. tout( kl ) ) then
      			fout( k ) = fout( k ) + xmat( kl )
      		end if
      	end do
      end do

      return
      end


! -- Pagebreak --
! Pg. 112:


      ! Note: While the formating code in this routine has been carefully reproduced the output
      ! still has some alignment issues.
      subroutine SUMMARY( outfil , number , maxno , maxkl, parts , nun ,
     +tis , ak , wdry , fmois , sigma , tign , tout , xmat , wo , diam )

      real*4 wdry( maxno) , fmois( maxno ) , sigma( maxno )
      real*4 tign( maxkl ) , tout( maxkl ) , wo( maxkl )
      real*4 xmat( maxkl ) , diam( maxkl )
      character*12 parts( maxno ) , outfil , name , none , nuname
      character*3 stat
      logical v

      loc( k , l ) = k * ( k + 1 ) / 2 + l

      none = 'no companion'
      nuname = outfil
      stat = 'NEW'

01    write( * , 05 )
05    format(' Full summary = 1 , recap only = 0    ',$ )
      read( * , * ) ido
      v = ( ido .EQ. 1 )
      if( ( ido .NE. 0 ) .AND. ( .NOT. v ) ) goto 01

10    open( nun , file=nuname , status=stat , form='FORMATTED' , err=99)
      write( nun , 15 ) ak , tis
15    format(5x'Area factor =' , f6.2 , '  End time , s =' , e12.4 )
20    format(/5x,a12,3x'Load , mois , diam =  '3e10.3/
     +5x'  Companion  fraction   Ignition   burnout   load rem    diam')
21    format(6x'Component   load @ 0   Ignition   burnout   load rem')
30    format(5x,a12,e9.3,e11.3,e10.3,e11.3,e9.3)
      if( .NOT. v ) write( nun , 21 )
35    do m = 1 , number
      	name = parts( m )
      	win = wdry( m )
      	fmi = fmois( m )
      	dim = 4. / sigma( m )
      	if( v ) write( nun , 20 ) name , win , fmi , dim
      	name = none
      	rem = 0.
      	ts = 1.e+31
      	tf = 0.
      	do n = 0 , m
      		mn = loc( m , n )
      		fr = xmat( mn )
      		ti = tign( mn )
      		ts = min( ti , ts )
      		to = tout( mn )
      		tf = max( to , tf )
      		wd = wo ( mn )
      		rem = rem + wd
      		di = diam( mn )
      		if( v ) write( nun , 30 ) name , fr , ti , to , wd , di
      		if( n .LT. m ) name = parts( ( n + 1 ) )
      	end do
      	name = parts( m )
      	if( .NOT. v) write( nun, 40 ) name , win , ts , tf , rem


! -- Pagebreak --
! Pg. 113:


      end do
40    format(5x,a12,e9.3,e11.3,e10.3,e11.3)
      return
99    write(*,100) nuname , stat
100   format('  Error on open:  'a12'  with status = 'a3/
     +     '  0 = abort now'/
     +     '  1 = overwrite existing file'/
     +     '  2 = enter new file name to create   ',$)
      read(*,*) in
      if( in .EQ. 0) stop' Abort'
      if( in .EQ. 1 ) then
      	stat = 'OLD'
      	goto 10
      end if
      if( in .NE. 2 ) goto 99
      write(*,*)' Enter new file name [max 12 characters]'
      read(*,101) nuname
101   format(a12)
      stat = 'NEW'
      open( nun , file=nuname , status=stat , form='FORMATTED' , err=99)
      goto 35
      ! JMR_NOTE: I believe this goto is a bug in the original code.  It will cause the first line,
      ! with the area factor and end time, to be omitted from the summary file.  I think the code
      ! should go to the write() on the line after 10 instead.
      
      end


! End of source code.
