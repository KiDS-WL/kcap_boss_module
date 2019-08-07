import os
import ctypes as ct
import numpy as np

def array_ctype(ndim, dtype=np.float64, flags="F"):
    return [np.ctypeslib.ndpointer(ndim=ndim, dtype=dtype, flags=flags)] + [ct.POINTER(ct.c_int)]*ndim

def array_arg(a):
    arr = a
    return (arr, *(ct.c_int(s) for s in arr.shape))

class BiasModule:
    libname = "lib/libbias_wrapper.so"
    module_name = "bias_module"

    def __init__(self):
        self.load_lib()

    def load_lib(self, path=None):
        if path is None:
            path = os.path.dirname(__file__)
        libpath = os.path.abspath(os.path.join(path, self.libname))
        self.lib = ct.CDLL(libpath)

    def get_function(self, name, c_bind=True):
        if c_bind:
            return getattr(self.lib, name)
        else:
            return getattr(self.lib, f"__{self.module_name}_MOD_{name}")

    def initialize_wedges(self, twopt_type, num_ell, num_points_use, num_bands_use, zm, om_fid, h0_fid, 
                                window,
                                verbose=True):
        f = self.get_function("initialize_wedges")
        f.restype = ct.c_void_p
        f.argtypes = [ct.POINTER(ct.c_int),        # twopt_type
                      ct.POINTER(ct.c_int),        # num_ell
                      ct.POINTER(ct.c_int),        # num_points_use
                      ct.POINTER(ct.c_int),        # num_bands_use
                      ct.POINTER(ct.c_double),     # zm
                      ct.POINTER(ct.c_double),     # om_fid
                      ct.POINTER(ct.c_double),     # h0_fid
                      *array_ctype(ndim=2, dtype=np.float64), # window
                      ct.POINTER(ct.c_int),        # verbose
                      ]      

        dataset_ptr = f(ct.c_int(twopt_type), ct.c_int(num_ell), ct.c_int(num_points_use), 
                        ct.c_int(num_bands_use),
                        ct.c_double(zm), ct.c_double(om_fid), ct.c_double(h0_fid),
                        *array_arg(np.asfortranarray(window)),
                        ct.c_int(verbose))
        self.num_points_use = num_points_use
        self.num_ell = num_ell
        return dataset_ptr

    def cleanup_wedges(self, config_ptr):
        f = self.get_function("cleanup_wedges")
        f.restype = None
        f.argtypes = [ct.c_void_p]
        f(config_ptr)

    def setup_cosmology(self, h, omdm, omb, omv, omk, omnuh2, nnu, w, wa,
                              use_growth, local_lag_g2, local_lag_g3, 
                              z, log_k, Pk, 
                              growth_z, sigma8_z, 
                              verbose=True):
        f = self.get_function("setup_cosmology")
        f.restype = None
        f.argtypes = [ct.POINTER(ct.c_double),     # h
                      ct.POINTER(ct.c_double),     # omdm
                      ct.POINTER(ct.c_double),     # omb
                      ct.POINTER(ct.c_double),     # omv
                      ct.POINTER(ct.c_double),     # omk
                      ct.POINTER(ct.c_double),     # omnuh2
                      ct.POINTER(ct.c_double),     # nnu
                      ct.POINTER(ct.c_double),     # w
                      ct.POINTER(ct.c_double),     # wa
                      ct.POINTER(ct.c_bool),       # use_growth
                      ct.POINTER(ct.c_bool),       # local_lag_g2
                      ct.POINTER(ct.c_bool),       # local_lag_g3
                      *array_ctype(ndim=1, dtype=np.float64), # Pk_z
                      *array_ctype(ndim=1, dtype=np.float64), # Pk_log_k
                      *array_ctype(ndim=2, dtype=np.float64), # Pk
                      *array_ctype(ndim=1, dtype=np.float64), # growth_z
                      *array_ctype(ndim=1, dtype=np.float64), # sigma8_z
                      ct.POINTER(ct.c_int),        # verbose
                      ]      

        f(ct.c_double(h), ct.c_double(omdm), ct.c_double(omb), ct.c_double(omv), ct.c_double(omk),
          ct.c_double(omnuh2), ct.c_double(nnu), ct.c_double(w), ct.c_double(wa),
          ct.c_bool(use_growth), ct.c_bool(local_lag_g2), ct.c_bool(local_lag_g3), 
          *array_arg(z),
          *array_arg(log_k),
          *array_arg(np.asfortranarray(Pk)),
          *array_arg(growth_z),
          *array_arg(sigma8_z),
          ct.c_int(verbose)
          )
        self.log_k = log_k

    def compute_wedges(self, config_ptr,
                             b1, b2, gamma2, gamma3, a_vir, gamma,
                             z_index,
                             H_z, DA_z, 
                             bands, 
                             verbose=True):
        f = self.get_function("compute_wedges")
        f.restype = None
        f.argtypes = [ct.c_void_p,
                      ct.POINTER(ct.c_double),     # b1
                      ct.POINTER(ct.c_double),     # b2
                      ct.POINTER(ct.c_double),     # gamma2
                      ct.POINTER(ct.c_double),     # gamma3
                      ct.POINTER(ct.c_double),     # a_vir
                      ct.POINTER(ct.c_double),     # gamma
                      ct.POINTER(ct.c_int),        # z_index
                      *array_ctype(ndim=1, dtype=np.float64), # H_z
                      *array_ctype(ndim=1, dtype=np.float64), # DA_z
                      *array_ctype(ndim=1, dtype=np.float64), # bands
                      *array_ctype(ndim=1, dtype=np.float64), # vtheo
                      *array_ctype(ndim=1, dtype=np.float64), # vtheo_convolved
                      *array_ctype(ndim=1, dtype=np.float64), # Pk_mm
                      *array_ctype(ndim=1, dtype=np.float64), # Pk_gm
                      *array_ctype(ndim=1, dtype=np.float64), # Pk_gg
                      ct.POINTER(ct.c_int),        # verbose
                      ]      
        vtheo = np.empty(self.num_ell*len(bands))
        vtheo_convolved = np.empty(self.num_ell*self.num_points_use)

        Pk_mm = np.zeros(len(self.log_k))
        Pk_gm = np.zeros(len(self.log_k))
        Pk_gg = np.zeros(len(self.log_k))

        f(config_ptr, 
          ct.c_double(b1), ct.c_double(b2), ct.c_double(gamma2), ct.c_double(gamma3), ct.c_double(a_vir), ct.c_double(gamma),
          ct.c_int(z_index), 
          *array_arg(np.asfortranarray(H_z)), *array_arg(np.asfortranarray(DA_z)),
          *array_arg(bands),
          *array_arg(vtheo),
          *array_arg(vtheo_convolved),
          *array_arg(Pk_mm),
          *array_arg(Pk_gm),
          *array_arg(Pk_gg),
          ct.c_int(verbose)
          )
        return vtheo, vtheo_convolved, Pk_mm, Pk_gm, Pk_gg


if __name__ == "__main__":
    h = 0.67956542968750000     
    omdm = 0.25731178759694234     
    omb = 4.8246962797075028e-002
    omv = 0.69444124960598264
    omk = 0.0
    omnuh2 = 6.4514389153979819e-4
    nnu = 3.046
    w = -1.0
    wa = 0.0

    b1 = 2.069093
    b2 = -0.08266389
    gamma2 = 1.0
    gamma3 = 1.049944
    a_vir = 3.414406
    gamma = 0.47

    derived_parameters = np.loadtxt("../benchmark/derived_params.txt")
    H_z = derived_parameters[12+2::4]
    DA_z = derived_parameters[12+3::4]
    print("H(z)", H_z)
    print("DA(z)", DA_z)

    twopt_type = 4
    num_ell = 3
    num_points_use = 28
    num_bands_use = 140
    z_index = 4
    zm = 0.60999999999999999     
    om_fid = 0.31000000000000000     
    h0_fid = 0.69999999999999996

    use_growth = False 
    local_lag_g2 = True 
    local_lag_g3 = False

    window = np.loadtxt("../benchmark/window.txt")
    bands = np.loadtxt("../benchmark/bands.txt")
    
    z = np.loadtxt("../benchmark/z_p_k.txt")
    log_k = np.loadtxt("../benchmark/log_k_h.txt")
    Pk = np.loadtxt("../benchmark/log_p_k.txt")

    print(f"n_z: {len(z)}, n_k: {len(log_k)}, shape(Pk): {Pk.shape}")
    print(z)
    print(log_k[:5])
    print(Pk[:,0])

    growth_z = np.loadtxt("../benchmark/growth.txt")
    sigma8_z = np.loadtxt("../benchmark/sigma8.txt")

    mod = BiasModule()
    config_ptr = mod.initialize_wedges(twopt_type, num_ell, num_points_use, num_bands_use, zm, om_fid, h0_fid, 
                                       window, verbose=True)

    mod.setup_cosmology(h, omdm, omb, omv, omk, omnuh2, nnu, w, wa,
                        use_growth, local_lag_g2, local_lag_g3,
                        z, log_k, Pk, 
                        growth_z, sigma8_z)

    vtheo, vtheo_convolved, Pk_mm, Pk_gm, Pk_gg = mod.compute_wedges(config_ptr,
                             b1, b2, gamma2, gamma3, a_vir, gamma,
                             z_index,
                             H_z, DA_z, 
                             bands, 
                             )

    mod.cleanup_wedges(config_ptr)

    np.savetxt("../output/vtheo_python.txt", vtheo)
    np.savetxt("../output/vtheo_convolved_python.txt", vtheo_convolved)
