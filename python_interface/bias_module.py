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

    def compute_wedges(self, h, omdm, omb, omv, omk, omnuh2, nnu, w, wa, 
                             H_z, DA_z, 
                             twopt_type, num_ell, num_points_use, num_bands_use, z_index, zm, om_fid, h0_fid, 
                             use_growth, local_lag_g2, local_lag_g3, 
                             window, 
                             bands, 
                             z, log_k, Pk, 
                             growth_z, sigma8_z, 
                             data_params,
                             verbose=True):
        f = self.get_function("compute_wedges")
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
                      *array_ctype(ndim=1, dtype=np.float64), # H_z
                      *array_ctype(ndim=1, dtype=np.float64), # DA_z
                      ct.POINTER(ct.c_int),        # twopt_type
                      ct.POINTER(ct.c_int),        # num_ell
                      ct.POINTER(ct.c_int),        # num_points_use
                      ct.POINTER(ct.c_int),        # num_bands_use
                      ct.POINTER(ct.c_int),        # z_index
                      ct.POINTER(ct.c_double),     # zm
                      ct.POINTER(ct.c_double),     # om_fid
                      ct.POINTER(ct.c_double),     # h0_fid
                      ct.POINTER(ct.c_bool),       # use_growth
                      ct.POINTER(ct.c_bool),       # local_lag_g2
                      ct.POINTER(ct.c_bool),       # local_lag_g3
                      *array_ctype(ndim=2, dtype=np.float64), # window
                      *array_ctype(ndim=1, dtype=np.float64), # bands
                      *array_ctype(ndim=1, dtype=np.float64), # Pk_z
                      *array_ctype(ndim=1, dtype=np.float64), # Pk_log_k
                      *array_ctype(ndim=2, dtype=np.float64), # Pk
                      *array_ctype(ndim=1, dtype=np.float64), # growth_z
                      *array_ctype(ndim=1, dtype=np.float64), # sigma8_z
                      *array_ctype(ndim=1, dtype=np.float64), # data_params
                      *array_ctype(ndim=1, dtype=np.float64), # vtheo
                      *array_ctype(ndim=1, dtype=np.float64), # vtheo_convolved
                      ct.POINTER(ct.c_int),        # verbose
                      ]      
        vtheo = np.empty(num_ell*len(bands))
        vtheo_convolved = np.empty(num_ell*num_points_use)

        f(ct.c_double(h), ct.c_double(omdm), ct.c_double(omb), ct.c_double(omv), ct.c_double(omk),
          ct.c_double(omnuh2), ct.c_double(nnu), ct.c_double(w), ct.c_double(wa),
          *array_arg(np.asfortranarray(H_z)), *array_arg(np.asfortranarray(DA_z)),
          ct.c_int(twopt_type), ct.c_int(num_ell), ct.c_int(num_points_use), ct.c_int(num_bands_use), ct.c_int(z_index), 
          ct.c_double(zm), ct.c_double(om_fid), ct.c_double(h0_fid),
          ct.c_bool(use_growth), ct.c_bool(local_lag_g2), ct.c_bool(local_lag_g3), 
          *array_arg(np.asfortranarray(window)),
          *array_arg(bands),
          *array_arg(z),
          *array_arg(log_k),
          *array_arg(np.asfortranarray(Pk)),
          *array_arg(growth_z),
          *array_arg(sigma8_z),
          *array_arg(data_params),
          *array_arg(vtheo),
          *array_arg(vtheo_convolved),
          ct.c_int(verbose)
          )
        return vtheo, vtheo_convolved


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

    data_params = np.loadtxt("../benchmark/data_params.txt")

    mod = BiasModule()
    vtheo, vthep_convolved = mod.compute_wedges(
                             h, omdm, omb, omv, omk, omnuh2, nnu, w, wa, 
                             H_z, DA_z, 
                             twopt_type, num_ell, num_points_use, num_bands_use, z_index, zm, om_fid, h0_fid, 
                             use_growth, local_lag_g2, local_lag_g3,
                             window, 
                             bands, 
                             z, log_k, Pk, 
                             growth_z, sigma8_z, 
                             data_params)

    np.savetxt("../output/vtheo_python.txt", vtheo)
    np.savetxt("../output/vtheo_convolved_python.txt", vthep_convolved)
