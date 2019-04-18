

module fiducial
  use settings,  only: mcp        ! real precision
  use constants, only: const_pi  ! pi
  use numrec
  implicit none
  !values of dm and H used to take into account the fiducial cosmology. 
  !Assumes that the fiducial cosmology corresponds to a simple flat LCDM model 
  !as the one used for all BOSS measurements.
  real(mcp)       :: save_om_m !used to pass Omega_m to the integrand of D_M(z)
  contains

    subroutine compute_scaling(om_m_fid, h0_fid, z, d_m, hub)!Ariel: added this subroutine
      use camb
      implicit none
      real(kind=mcp), intent(in)  :: om_m_fid, h0_fid, z
      real(kind=mcp), intent(out) :: d_m, hub
      hub  = e(om_m_fid, z)!*h0_fid  !H(z)/100.
      d_m  = dm(om_m_fid, z)!/h0_fid !D_M(z)
    end subroutine compute_scaling
  
    function e(om_m, z)    !E(z)=H(z)/H_0
      implicit none
      real(kind=mcp)   :: om_m, z
      real(kind=mcp)   :: e2, e, zplus, a
      zplus  = 1._mcp + z
      a      = 1._mcp / zplus
      e2 = om_m/a**3 + (1._mcp- om_m)
      e = sqrt(e2)
    end function e

    function einv(z)
      implicit none
      real(kind=mcp) :: einv, z
      einv = 1._mcp/e(save_om_m, z)
    end function einv
  
    function dm(om_m, z) !for a flat cosmology dm = r(z)
      implicit none
      real(kind=mcp)   :: dm, om_m, z, ans
      save_om_m = om_m
      call qromb(einv, 0._mcp, z, ans)
      dm = ans*2997.92458_mcp
    end function dm

end module fiducial


