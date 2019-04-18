!Module to handle all two-point LSS datasets, their points and window functions
!and routines for computing the likelihood

!================================================================================
!: module containing a single type used to read and store all two-point LSS datasets
  module twopt_data_def
    use settings
    use CosmologyTypes
    use CosmoTheory
    use Calculator_Cosmology 
    use MatrixUtils
    !use Likelihood_Cosmology
    implicit none

    !: container for all the sizes and indices needed when reading the
    !dataset
    type :: sizes
      ! number of multipoles or wedges
      integer :: num_ell
      ! points: where the LSS statistic is measured
      ! bands: where to evaluate the model before convonving with the window matrix
      integer :: num_points_use ! total number of points used (ie. max-min+1)
      integer :: num_bands_use ! total number of kbands used (ie. max-min+1)
      !: numbers of points and kbands in the input files
      integer :: num_points_full 
      integer :: num_bands_full
      !: maximum and minimum indices of the measurements to use
      integer :: max_points_use, min_points_use 
      !: maximum and minimum indices where to evaluate the model before
      !convolving with the window matrix
      integer :: max_bands_use, min_bands_use 
    end type sizes

    type :: twopt_dataset
      character(LEN=80)                      :: name  !dataset name
      !which type of LSS measurements: pk_ell (1), pk_wed (2), xi_ell (3) or xi_wed (4)
      integer                                :: twopt_type      
      type(sizes)                            :: nsize         !storage of the sizes and indeced to use
      real(mcp), allocatable, dimension(:)   :: measur        !P(k), xi0, xi_ell or xi_wed
      real(mcp), allocatable, dimension(:)   :: points        !scales corresponding to the measurements
      real(mcp), allocatable, dimension(:,:) :: invcov        !inverse covariance
      real(mcp), allocatable, dimension(:)   :: bands         !scales where the model is evaluated before convolving
      real(mcp), allocatable, dimension(:,:) :: window        !window matrix 
      real(mcp), allocatable, dimension(:)   :: zerowindowfxn !W(k=0, kj), only used for P(k)
      !window function evaluated in k and in k=0 (the scalar)
      real(mcp), allocatable, dimension(:)   :: zerowindowfxnsubtractdat !only used for P(k)
      real(mcp)                              :: zerowindowfxnsubdatnorm  !only used for P(k)
      real(mcp)                              :: zm        !mean redshift of the sample
      real(mcp)                              :: om_fid, h0_fid  
      contains
        !methods
        procedure          :: read_data      ! public procedure to read the dataset
        procedure, private :: read_sizes     ! get the indices and sizes
        procedure, private :: read_measurements   ! read two-point LSS measurement
        procedure, private :: read_invert_covmat  ! read and invert the covariance matrix
        procedure, private :: read_window         ! read the files of the window function
        procedure, private :: read_fiducial       ! read and set the fiducial cosmological parameters
        !convolve input power spectrum with the window matrix
        !procedure, private :: average_xi         !bin average of the correlation function
        !procedure, private :: convolve_pk        !window matrix convolution of the model P(k)
        !procedure, pointer :: convolve => null() !pointer to convolve_pk or average_xi
        procedure          :: convolve           !convolve pk or average xi
    end type twopt_dataset

    contains !: method definitions

    !call the methods that get all the
    !important data from the 'dataset' file
    subroutine read_data(mset, Ini)
      implicit none
      class(twopt_dataset)                           :: mset
      class(TSettingIni)                             :: Ini
      character(len=20)                              :: data_type

      mset%name = Ini%Read_String('name')     !dataset name
      if(feedback > 0) write(*,*) 'Reading dataset: '//trim(mset%name)
      !get dataset type
      data_type = Ini%Read_String('data_type')     !dataset name
      if(trim(data_type)=='pk_ell')then
        mset%twopt_type=1
        !mset%convolve => mset%convolve_pk
      elseif(trim(data_type)=='pk_wed')then
        mset%twopt_type=2
        !mset%convolve => mset%average_xi
      elseif(trim(data_type)=='xi_ell')then
        mset%twopt_type=3
        !mset%convolve => mset%average_xi
      elseif(trim(data_type)=='xi_wed')then
        mset%twopt_type=4
        !mset%convolve => mset%average_xi
      else
        write(*,*)'incorrect LSS data type in '//trim(mset%name)
        stop
      end if
      !get dataset name

      !get the indices and sizes
      call mset%read_sizes(Ini)
      !read the measurements
      call mset%read_measurements(Ini)
      !read and invert the covariance matrix
      call mset%read_invert_covmat(Ini)
      !read the files of the window
      call mset%read_window(Ini)
      !fiducial cosmology and mean redshift
      call mset%read_fiducial(Ini)
    end subroutine read_data

    !get the indices and sizes to read from the dataset files
    subroutine read_sizes(mset, Ini)
      implicit none
      class(twopt_dataset)         :: mset  !object reference
      class(TSettingIni)         :: Ini
      !the total number of multipoles or wedges  must be given
      mset%nsize%num_ell         = Ini%Read_Int('num_ell',0)
      if (mset%nsize%num_ell == 0) write(*,*) ' ERROR: parameter num_ell not set'
      !the total number of points and bands must be given
      mset%nsize%num_points_full = Ini%Read_Int('num_points_full',0)
      if (mset%nsize%num_points_full == 0) write(*,*) ' ERROR: parameter num_points_full not set'
      mset%nsize%num_bands_full  = Ini%Read_Int('num_bands_full',0)
      if (mset%nsize%num_bands_full == 0) write(*,*) ' ERROR: parameter num_bands_full not set'
      !if any of these is not present,  defaults to 1 for 'min_*' and 'num_*_full' for 'max_*'
      mset%nsize%min_points_use  = Ini%Read_Int('min_points_use', 1)
      mset%nsize%min_bands_use   = Ini%Read_Int('min_bands_use', 1)
      mset%nsize%max_points_use  = Ini%Read_Int('max_points_use', mset%nsize%num_points_full)
      mset%nsize%max_bands_use   = Ini%Read_Int('max_bands_use', mset%nsize%num_bands_full)
      !actual numbers of bins used
      mset%nsize%num_points_use  = mset%nsize%max_points_use - mset%nsize%min_points_use +1
      mset%nsize%num_bands_use   = mset%nsize%max_bands_use - mset%nsize%min_bands_use +1
    end subroutine read_sizes

    !read the measured power spectrum
    subroutine read_measurements(mset, Ini)
      implicit none
      class(twopt_dataset)                 :: mset  !object reference
      class(TSettingIni)                   :: Ini
      character(180)                       :: measurements_file !: LSS file name
      character(180)                       :: readchar   !: store each line of the file
      real(mcp)                            :: dummy      !: dummy variable 
      integer                              :: i, j, k    !: loop variables
      integer                              :: ndat, nell   !: 
      integer                              :: statio     !: status variable when reading the file
      Type(TTextFile)                      :: F

      !read the measurement file
      measurements_file  = Ini%ReadFileName('measurements_file')
      if(feedback > 1) write(*,*) 'Reading file:'//trim(measurements_file)
      !allocate the arrays containing measurements
      allocate(mset%points(mset%nsize%num_points_use)) 
      !the total number of data points is given by num_ell*num_points_use 
      allocate(mset%measur(mset%nsize%num_ell*mset%nsize%num_points_use))
      !allocate auxiliary arrays for reading data sets
      ndat = mset%nsize%num_points_use
      nell = mset%nsize%num_ell

      mset%measur = 0._mcp
      !open, read and close the file
      call F%Open(measurements_file)

      i = 0 !file lines counter
      j = 0 !file lines counter
      do  !read line by line
        read(F%unit, '(a)', iostat=statio) readchar  !: read a line
        if(statio < 0) exit  !end of file reached
        if(statio > 0) then  !error occurred
          write (*,*) 'Error: ', statio, ' occurred while reading file ',trim(measurements_file)
          stop
        end if
        if(scan(readchar, '#') /= 0) cycle  !header, skip
        i = i + 1   !the line is not a comment, increase counter
        if(i < mset%nsize%min_points_use) cycle !first points not used
        j = j + 1   !the line must be read and stored

        read(readchar, *) mset%points(j), (mset%measur(j+(k-1)*ndat), dummy, k=1, nell)

        if(j == mset%nsize%num_points_use) exit !all required data points have been read and stored
      end do  !end of while loop
      call F%Close()

      if(feedback > 1) write(*,*) j, ' bins of twopt LSS data read'
    end subroutine read_measurements

    !read and invert the covariance matrix
    subroutine read_invert_covmat(mset, Ini)
      implicit none
      class(twopt_dataset)                        :: mset     !object reference
      class(TSettingIni)                     :: Ini
      character (LEN=180)                    :: cov_file !covariance matrix file name
      integer                                :: n_mocks  !number of mocks used to estimate the covariance matrix
      real(mcp)                              :: cor_fact !correction of the inverse covariance matrix
      real(mcp), allocatable, dimension(:)   :: cov_temp !temporary storage for the covariance
      character(10)                          :: readchar !used to check for header lines
      integer                                :: statio   !status variable when reading the file
      integer                                :: nhead, i, j, k  !loop variables
      integer                                :: nell, ndat, nfull, nmin, nmax
      integer                                :: istart, iend  
      logical                                :: use_line
      Type(TTextFile)                        :: F


      !read the file name and number of mocks
      cov_file = Ini%ReadFileName('cov_file')
      n_mocks  = Ini%Read_Int('n_mocks', 2048)
      if(feedback > 1) write(*,*) 'Reading file:'//trim(cov_file)
      !allocate the covariance and the inverse
      ndat  = mset%nsize%num_points_use
      nmin  = mset%nsize%min_points_use
      nmax  = mset%nsize%max_points_use
      nfull = mset%nsize%num_points_full
      nell  = mset%nsize%num_ell

      allocate(cov_temp(nell*nfull))
      allocate(mset%invcov(nell*ndat, nell*ndat))

      !read the covariance matrix
      call F%Open(cov_file)
      nhead = 0
      do  !read line by line
        read(F%unit, '(a)', iostat=statio) readchar  !: read a line
        if(statio < 0) exit  !end of file reached
        if(statio > 0) then  !error occurred
          write (*,*) 'Error: ', statio, ' occurred where reading file ',&
            trim(cov_file)
          stop
        end if
        if(scan(readchar, '#') /= 0) then
          nhead = nhead + 1
          cycle  !header, skip line
        end if
        exit
      end do
      write(*,*)'Header contains ',nhead, ' lines'
      rewind(F%unit)

      do i = 1,nhead
        read(F%unit, '(a)', iostat=statio) readchar  !: read a line
      end do

      i = 0
      j = 0
      do  !read line by line
        read(F%unit, *, iostat=statio) cov_temp  !: read a line
        if(statio < 0) exit  !end of file reached
        if(statio > 0) then  !error occurred
          write (*,*) 'Error: ', statio, ' occurred where reading file ',&
            trim(cov_file)
          stop
        end if
        i = i + 1   !the line is not a comment, increase counter
        !read row of full covariance matrix
        use_line = .false.
        do k = 1, nell
          use_line = use_line .or. ((i >= (nfull*(k-1) + nmin)) .and. (i <= (nfull*(k-1)+nmax)))
        end do
        if(use_line)then
          j = j + 1
          do k = 1, nell
            istart = ndat*(k-1) + 1
            iend   = ndat*k
            mset%invcov(j,istart:iend) = cov_temp(nfull*(k-1) + nmin : nfull*(k-1) + nmax)
          end do
        end if
      end do  
      if(i /= nfull*nell)then
        write (*,*) 'Error: covariance matrix file ended before C was fully read'
        stop
      end if

      call F%Close()

      deallocate(cov_temp)
      !invert the covariance matrix
      call Matrix_Inverse(mset%invcov)

      !Hartlap et al. (2007): the inverse of the maximum-likelihood
      !estimator of the covariance is biased.
      !inverse assumed unbiased if n_mocks <= 0
      if(n_mocks > 0) then
        !: if n-p-2<0 it is not possible to unbias the inverse
        if(n_mocks < ndat + 2) then  
          write(*,*) 'Insuficient number of mock catalogues.'
          write(*,*) 'The covariance matrix is singular'
          stop
        else
          cor_fact = (real(nell*ndat,kind=mcp) + 1._mcp)/(real(n_mocks,kind=mcp) - 1._mcp)
          mset%invcov = mset%invcov*(1._mcp - cor_fact)
        end if
      end if
      if(feedback > 1) write(*,*) 'Covariance matrix read and inverted'
    end subroutine read_invert_covmat

    !read window function
    subroutine read_window(mset, Ini)
      implicit none
      integer                                :: i, j
      integer                                :: istart, iend, jstart, jend
      integer                                :: iistart, iiend, jjstart, jjend
      class(twopt_dataset)                   :: mset  !object reference
      class(TSettingIni)                     :: Ini
      real(mcp), allocatable, dimension(:)   :: temp_v !: temporary vectore storage
      real(mcp), allocatable, dimension(:,:) :: temp_m !: temporary matrix storage
      character(LEN=180)                     :: bands_file,& !: bands file
        zerowindowfxn_file, zerowindowfxnsubtractdat_file, & !: integral constraint
        windows_file !: windo matrix

      !allocate and read kbands_file
      allocate(mset%bands(mset%nsize%num_bands_use)) !kj
      allocate(temp_v(mset%nsize%num_bands_full))    !temporary storage

      bands_file  = Ini%ReadFileName('bands_file')      !get file name
      if(feedback > 1) write(*,*) 'Reading file:'//trim(bands_file)
      call File%ReadTextVector(bands_file, temp_v, mset%nsize%num_bands_full)
      mset%bands=&  !copy only the interesting part of the vector
        temp_v(mset%nsize%min_bands_use:mset%nsize%max_bands_use) 

      if(mset%twopt_type == 1 .or. mset%twopt_type == 2)then
        !allocate and read the window matrix
        allocate(mset%window(mset%nsize%num_ell*mset%nsize%num_points_use, mset%nsize%num_ell*mset%nsize%num_bands_use))
        allocate(temp_m(mset%nsize%num_ell*mset%nsize%num_points_full, mset%nsize%num_ell*mset%nsize%num_bands_full))
        windows_file  = Ini%ReadFileName('windows_file')
        call File%ReadTextMatrix(windows_file, temp_m, mset%nsize%num_ell*mset%nsize%num_points_full, &
                                & mset%nsize%num_ell*mset%nsize%num_bands_full)
        do i = 1, mset%nsize%num_ell
          istart = mset%nsize%num_points_use*(i-1) + 1
          iend   = mset%nsize%num_points_use*i
          iistart = mset%nsize%num_points_full*(i-1) + mset%nsize%min_points_use
          iiend   = mset%nsize%num_points_full*(i-1) + mset%nsize%num_points_use
          do j = 1, mset%nsize%num_ell
            jstart = mset%nsize%num_bands_use*(j-1) + 1
            jend   = mset%nsize%num_bands_use*j
            jjstart = mset%nsize%num_bands_full*(j-1) + mset%nsize%min_bands_use
            jjend   = mset%nsize%num_bands_full*(j-1) + mset%nsize%num_bands_use
            mset%window(istart:iend,jstart:jend) = temp_m(iistart:iiend, jjstart:jjend) 
          end do
        end do
        deallocate(temp_m)
      else if(mset%twopt_type == 3 .or. mset%twopt_type == 4)then 
        !allocate and read the window matrix
        allocate(mset%window(mset%nsize%num_points_use, mset%nsize%num_bands_use))
        allocate(temp_m(mset%nsize%num_points_full, mset%nsize%num_bands_full))
        windows_file  = Ini%ReadFileName('windows_file')
        call File%ReadTextMatrix(windows_file, temp_m, mset%nsize%num_points_full, mset%nsize%num_bands_full)
        mset%window = temp_m(mset%nsize%min_points_use:mset%nsize%max_points_use,&
        mset%nsize%min_bands_use:mset%nsize%max_bands_use)
        deallocate(temp_m)
      end if

      if(feedback > 1) write(*,*) 'Window function read'
    end subroutine read_window

    !read fiducial cosmology, assuming a flat LCDM
    subroutine read_fiducial(mset, Ini)
      implicit none
      class(twopt_dataset)                   :: mset  !object reference
      class(TSettingIni)                     :: Ini
      real(kind=mcp)                         :: om, ol, hh
      mset%om_fid  = Ini%Read_Double('omega_fid', 0.31_mcp)  
      mset%h0_fid  = Ini%Read_Double('h0_fid', 0.7_mcp)  
      mset%zm = Ini%Read_Double('mean_redshift', 0.57_mcp)  
      if(feedback > 1)then
        write(*,*)'Omega_fid = ', mset%om_fid
        write(*,*)'h0_fid = ', mset%h0_fid
        write(*,*)'mean redshift = ', mset%zm
      end if
    end subroutine read_fiducial

    !Convolves the input model with the window functions.
    subroutine convolve(mset, model_in, model_out)
      implicit none
      class(twopt_dataset)           :: mset         !object reference
      real(kind=mcp), intent(in)     :: model_in(:)  !input P(k) to convolve
      real(kind=mcp), intent(out)    :: model_out(:) !output convolved P(k)
      real(kind=mcp)                 :: normW0       !normalisation of the integral constrain
      integer                        :: i, istart, iend, jstart, jend
      if(mset%twopt_type <= 2)then
        !for Fourier space measurements each wedge has its own window
        !write(*,*)'convolution with window function not implemented correctly!!!'
        model_out = matmul(mset%window,model_in)
      else
        !in configuration space al "windows" are the same and correspond to the distance binning
        do i = 1, mset%nsize%num_ell 
          istart = mset%nsize%num_points_use*(i-1) + 1
          iend   = mset%nsize%num_points_use*i
          jstart = mset%nsize%num_bands_use*(i-1) + 1
          jend   = mset%nsize%num_bands_use*i
          model_out(istart:iend) = matmul(mset%window,model_in(jstart:jend))
        end do
      end if
    end subroutine convolve

  end module twopt_data_def

