import os
import numpy as np
import scipy.interpolate

from cosmosis.datablock import option_section, names

from cosmosis.datablock.cosmosis_py import errors

import bias_module

def setup(options):
    config = {}
    
    config["module"] = bias_module.BiasModule()
    
    here = os.path.dirname(__file__)
    window_file = options[option_section, "window_file"]
    bands_file = options[option_section, "bands_file"]
    config["window"] = np.loadtxt(window_file)
    config["bands"] = np.loadtxt(bands_file)

    config["twopt_type"] = 4
    config["num_ell"] = 3
    config["num_points_use"] = 28
    config["num_bands_use"] = 140
    config["z_index"] = 4
    config["zm"] = 0.61
    config["om_fid"] = 0.31
    config["h0_fid"] = 0.7

    config["use_growth"] = False
    config["local_lag_g2"] = True
    config["local_lag_g3"] = False

    config["derived_parameters"] = np.loadtxt(os.path.join(here, "../output/derived_params.txt"))
    config["data_parameters"] = np.loadtxt(os.path.join(here, "../output/data_params.txt"))

    config["z_pk"] = np.loadtxt(os.path.join(here, "../output/z_p_k.txt"))
    config["log_k_pk"] = np.loadtxt(os.path.join(here, "../output/log_k_h.txt"))
    config["log_pk"] = np.loadtxt(os.path.join(here, "../output/log_p_k.txt"))

    return config

def execute(block, config):
    h = block[names.cosmological_parameters, "h0"]
    omdm = block[names.cosmological_parameters, "omega_c"]
    omb = block[names.cosmological_parameters, "omega_b"]
    omv = block[names.cosmological_parameters, "omega_lambda"]
    omk = block[names.cosmological_parameters, "omega_k"]
    omnuh2 = block[names.cosmological_parameters, "omnuh2"]
    nnu = block.get_int(names.cosmological_parameters, "massive_nu", default=0) \
         + block.get_double(names.cosmological_parameters, "massless_nu", default=3.046)
    w = block[names.cosmological_parameters, "w"]
    wa = block[names.cosmological_parameters, "wa"]

    z = block[names.distances, "z"]
    z_index = np.argmin(np.abs(z-config["zm"]))
    H_z = block[names.distances, "h"]*2.99792458e8/1e3
    DA_z = block[names.distances, "d_a"]

    if z_index >= (200-13)//4:
        # z_index too large to fit H_z and DA_z into derived_parameters array
        # so we take a range around z[z_index]
        s = slice(max(0, z_index-4), min(len(z), z_index+4))
        z = z[s]
        H_z = H_z[s]
        DA_z = DA_z[s]
        z_index = np.argmin(np.abs(z-config["zm"]))

    print(f"z_index: {z_index} (z={z[z_index]})")
    print(f"H(z_index) = {H_z[z_index]}, DA(z_index) = {DA_z[z_index]}")


    log_Pk = np.log(block[names.matter_power_lin, "p_k"])
    log_k_h = np.log(block[names.matter_power_lin, "k_h"])
    z = block[names.matter_power_lin, "z"]

    z_growth = block["growth", "z"]
    if not np.allclose(z, z_growth):
        raise ValueError("Redshifts of power spectrum and growth do not match.")

    sigma8 = block["growth", "SIGMA_8"]
    sigma2_vdelta_8 = block["growth", "SIGMA2_VDELTA_8"]
    growth = sigma2_vdelta_8/sigma8

    # log_Pk_intp = scipy.interpolate.RectBivariateSpline(z, log_k_h, log_Pk)
    # log_k_h = config["log_k_pk"]
    # z = config["z_pk"]
    # log_Pk = log_Pk_intp(z, log_k_h, grid=True)

    # log_Pk = config["log_pk"].T
    # log_k_h = config["log_k_pk"]
    # z = config["z_pk"]

    # sigma8 = scipy.interpolate.InterpolatedUnivariateSpline(z_growth, sigma8)(z)
    # growth = scipy.interpolate.InterpolatedUnivariateSpline(z_growth, growth)(z)
    # z_growth = z

    vtheo, vtheo_convolved = config["module"].compute_wedges(
                             h, omdm, omb, omv, omk, omnuh2, nnu, w, wa, 
                             H_z, DA_z,
                             config["twopt_type"], config["num_ell"], 
                             config["num_points_use"], config["num_bands_use"], 
                             z_index+2, config["zm"], config["om_fid"], config["h0_fid"], 
                             config["use_growth"], config["local_lag_g2"], config["local_lag_g3"],
                             config["window"], 
                             config["bands"], 
                             z, log_k_h, log_Pk.T,
                             growth, sigma8,
                             config["data_parameters"],
                             verbose=True)

    n = len(config["bands"])
    vtheo = np.array([vtheo[i*n:(i+1)*n] for i in range(config["num_ell"])])
    n = config["num_points_use"]
    vtheo_convolved = np.array([vtheo_convolved[i*n:(i+1)*n] for i in range(config["num_ell"])])

    block["wedges", "vtheo"] = vtheo
    block["wedges", "bands"] = config["bands"]
    block["wedges", "vtheo_convolved"] = vtheo_convolved

    return 0

def cleanup(config):
    pass
