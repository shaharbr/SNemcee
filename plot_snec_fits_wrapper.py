import os
from matplotlib import pyplot as plt
import plot_snec_fits
# import plot_snec_fits_martinez_and_data
import sys, getopt



def lum_wCSM_vs_woCSM(SN_name, results_dir, output_dir):
    fig, axs = plt.subplots(1, 2, sharey='row', figsize=(20, 12))

    lum_csmTrue_name = SN_name + '_lum_csm-with_normalizedFalse_TreshLumFalse'
    lum_csmTrue_path = os.path.join(results_dir, lum_csmTrue_name)

    lum_csmFalse_name = SN_name + '_lum_csm-without_normalizedFalse_TreshLumFalse'
    lum_csmFalse_path = os.path.join(results_dir, lum_csmFalse_name)

    plot_snec_fits.plot_result_fit(lum_csmTrue_path, 'lum', axs[0])
    plot_snec_fits.plot_result_fit(lum_csmFalse_path, 'lum', axs[1])

    axs[0].set_ylabel('Log bolometric luminosity (erg/s)', fontsize=14)
    axs[0].set_xlabel('Rest-frame days from discovery', fontsize=14)
    axs[0].set_xlabel('Rest-frame days from discovery', fontsize=14)

    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, SN_name+'_lum_csm_comparison.png'))

    plot_snec_fits.overlay_corner_plot([lum_csmFalse_path, lum_csmTrue_path], output_dir,
                                       ['without CSM', 'with CSM'], SN_name + '_lum_csm_comparison')

    return fig


def lum_vs_lum_veloc_vs_lum_veloc_normalized(SN_name, results_dir, output_dir):
    fig, axs = plt.subplots(2, 3, sharey='row', figsize=(20, 12))

    lum_name = SN_name + '_lum_csm-with_normalizedFalse_TreshLumFalse'
    lum_path = os.path.join(results_dir, lum_name)
    lum_veloc_name = SN_name + '_lum-veloc_csm-with_normalizedFalse_TreshLumFalse_TreshVelocTrue'
    lum_veloc_path = os.path.join(results_dir, lum_veloc_name)
    lum_veloc_normalized_name = SN_name + '_lum-veloc_csm-with_normalizedTrue_TreshLumFalse_TreshVelocTrue'
    lum_veloc_normalized_path = os.path.join(results_dir, lum_veloc_normalized_name)

    plot_snec_fits.plot_result_fit(lum_path, 'lum', axs[0, 0])
    plot_snec_fits.plot_result_fit(lum_veloc_path, 'lum', axs[0, 1])
    plot_snec_fits.plot_result_fit(lum_veloc_normalized_path, 'lum', axs[0, 2])
    plot_snec_fits.plot_result_fit(lum_path, 'veloc', axs[1, 0])
    plot_snec_fits.plot_result_fit(lum_veloc_path, 'veloc', axs[1, 1])
    plot_snec_fits.plot_result_fit(lum_veloc_normalized_path, 'veloc', axs[1, 2])

    axs[0, 0].set_ylabel('Log bolometric luminosity (erg/s)', fontsize=14)
    axs[1, 0].set_ylabel('Expansion velocity (km/s)', fontsize=14)
    axs[1, 0].set_xlabel('Rest-frame days from discovery', fontsize=14)
    axs[1, 1].set_xlabel('Rest-frame days from discovery', fontsize=14)
    axs[1, 2].set_xlabel('Rest-frame days from discovery', fontsize=14)

    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, SN_name + '_lum_veloc_comparison.png'))

    plot_snec_fits.overlay_corner_plot([lum_path, lum_veloc_path, lum_veloc_normalized_path], output_dir,
                                       ['lum', 'lum+veloc, not normalized', 'lum+veloc, normalized'], SN_name + '_lum_veloc_comparison')
    return fig


def lum_veloc_vs_mag_veloc(SN_name, results_dir, output_dir, LumTthresh=False):
    fig, axs = plt.subplots(3, 2, sharey='row', figsize=(20, 12))

    if LumTthresh:
        lum_veloc_name = SN_name + '_lum-veloc_csm-with_normalizedFalse_TreshLumTrue_TreshVelocTrue'
    else:
        lum_veloc_name = SN_name + '_lum-veloc_csm-with_normalizedFalse_TreshLumFalse_TreshVelocTrue'
    lum_veloc_path = os.path.join(results_dir, lum_veloc_name)

    mag_veloc_name = SN_name + '_mag-veloc_csm-with_normalizedFalse_TreshMagTrue_TreshVelocTrue'
    mag_veloc_path = os.path.join(results_dir, mag_veloc_name)

    plot_snec_fits.plot_result_fit(lum_veloc_path, 'lum', axs[0, 0])
    plot_snec_fits.plot_result_fit(mag_veloc_path, 'lum', axs[0, 1])
    plot_snec_fits.plot_result_fit(lum_veloc_path, 'mag', axs[1, 0])
    plot_snec_fits.plot_result_fit(mag_veloc_path, 'mag', axs[1, 1])
    plot_snec_fits.plot_result_fit(lum_veloc_path, 'veloc', axs[2, 0])
    plot_snec_fits.plot_result_fit(mag_veloc_path, 'veloc', axs[2, 1])

    axs[0, 0].set_ylabel('Log bolometric luminosity (erg/s)', fontsize=14)
    axs[1, 0].set_ylabel('Absolute Magnitude', fontsize=14)
    axs[2, 0].set_ylabel('Expansion velocity (km/s)', fontsize=14)
    axs[2, 0].set_xlabel('Rest-frame days from discovery', fontsize=14)
    axs[2, 1].set_xlabel('Rest-frame days from discovery', fontsize=14)

    plt.tight_layout()
    fig.savefig(os.path.join(output_dir, SN_name + '_lum-veloc_mag-veloc_comparison_LumThresh'+str(LumTthresh)+'.png'))

    plot_snec_fits.overlay_corner_plot([lum_veloc_path, mag_veloc_path], output_dir,
                                       ['lum+veloc', 'mag+veloc'],
                                       SN_name + '_lum-veloc_mag-veloc_comparison_LumThresh'+str(LumTthresh))
    return fig


# TODO fix this one, after fiing run_mcmc_forall_twostep
def lum_veloc_onestep_vs_twostep(SN_name, results_dir, output_dir):
    fig, axs = plt.subplots(2, 3, sharey='row', figsize=(20, 12))

    onestep_name = SN_name + '_lum-veloc_csm-with_' + SN_name
    onestep_path = os.path.join(results_dir, onestep_name)
    twostep_priorNone_name = SN_name + '_lum-veloc_twostep_priorsNone_' + SN_name
    twostep_priorNone_path = os.path.join(results_dir, twostep_priorNone_name)
    twostep_priorTrue_name = SN_name + '_lum-veloc_twostep_priorsTrue_' + SN_name
    twostep_priorTrue_path = os.path.join(results_dir, twostep_priorTrue_name)


    plot_snec_fits.plot_result_fit(onestep_path, 'lum', axs[0, 0])
    plot_snec_fits.plot_result_fit(twostep_priorNone_path, 'lum', axs[0, 1])
    plot_snec_fits.plot_result_fit(twostep_priorTrue_path, 'lum', axs[0, 2])
    plot_snec_fits.plot_result_fit(onestep_path, 'veloc', axs[1, 0])
    plot_snec_fits.plot_result_fit(twostep_priorNone_path, 'veloc', axs[1, 1])
    plot_snec_fits.plot_result_fit(twostep_priorTrue_path, 'veloc', axs[1, 2])

    axs[0, 0].set_ylabel('Log bolometric luminosity (erg/s)', fontsize=14)
    axs[1, 0].set_ylabel('Expansion velocity (km/s)', fontsize=14)
    axs[1, 0].set_xlabel('Rest-frame days from discovery', fontsize=14)
    axs[1, 1].set_xlabel('Rest-frame days from discovery', fontsize=14)
    axs[1, 2].set_xlabel('Rest-frame days from discovery', fontsize=14)

    plt.tight_layout()

    fig.savefig(os.path.join(output_dir, SN_name + '_onestep_twostep_comparison.png'))

    plot_snec_fits.overlay_corner_plot([onestep_path, twostep_priorNone_path, twostep_priorTrue_path],
                                       output_dir,
                                       ['one step', 'two step, no priors carryover', 'two step, with priors carryover'],
                                       SN_name + '_onestep_twostep_comparison')
    return fig


def main(argv):
    SN_name = ""
    n_steps = 0
    type_fig = ""
    output_dir = ""
    arg_help = "{0} -S <SN name> -s <number of steps> " \
               "-t <type of figure [csm/lum-mag-veloc/lum-mag-veloc-Tthresh" \
               "/lum-veloc-normalized/lum-veloc-twostep]>" \
               " -o <output directory>".format(argv[0])
    try:
        opts, args = getopt.getopt(argv[1:], "hS:s:t:o:", ["help", "SN=",
                                                         "steps=", "type_figure=", "output_dir="])
    except:
        print(arg_help)
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print(arg_help)  # print the help message
            sys.exit(2)
        elif opt in ("-S", "--SN_name"):
            SN_name = arg
        elif opt in ("-s", "--steps"):
            n_steps = int(arg)
        elif opt in ("-t", "--type_figure"):
            type_fig = arg
        elif opt in ("-o", "--output_dir"):
            output_dir = arg

    print('SN name:', SN_name)
    print('steps:', n_steps)
    print('type of figure:', type_fig)
    print('output directory:', output_dir)
    res_dir = os.path.join('mcmc_results', output_dir, str(n_steps) + 'step')
    step_dir = os.path.join('figures', output_dir, str(n_steps) + 'step')
    if not os.path.exists(step_dir):
        if not os.path.exists(os.path.join('figures', output_dir)):
            os.mkdir(os.path.join('figures', output_dir))
        os.mkdir(step_dir)

    if type_fig == 'csm':
        lum_wCSM_vs_woCSM(SN_name, res_dir, step_dir)
    elif type_fig == 'lum-mag-veloc':
        lum_veloc_vs_mag_veloc(SN_name,res_dir, step_dir)
    elif type_fig == 'lum-mag-veloc-Tthresh':
        lum_veloc_vs_mag_veloc(SN_name,res_dir, step_dir, LumTthresh=True)
    elif type_fig == 'lum-veloc-normalized':
        lum_vs_lum_veloc_vs_lum_veloc_normalized(SN_name,res_dir, step_dir)
    elif type_fig == 'lum-veloc-twostep':
        lum_veloc_onestep_vs_twostep(SN_name,res_dir, step_dir)


if __name__ == "__main__":
    main(sys.argv)

# date = '2022_03_23'
# n_steps = 2000
# sn_names = ['SN2004a', 'SN2004et_d5.9', 'SN2012aw', 'SN2012ec', 'SN2017eaw','SN2018aoq', 'SN2005cs', 'SN2008bk']



# if not os.path.exists(os.path.join('mcmc_results', date)):
#     os.mkdir(os.path.join('mcmc_results', date))
# if not os.path.exists(os.path.join('mcmc_results', date, str(n_steps) + 'step')):
#     os.mkdir(os.path.join('mcmc_results', date, str(n_steps) + 'step'))
# if not res_dir:
#     os.mkdir(res_dir)
#
# if not os.path.exists(output_dir):
#     os.mkdir(output_dir)







# def lum_veloc_vs_martinez(SN_name, results_dir, output_dir):
#     fig, axs = plt.subplots(2, 3, sharey='row', figsize=(20, 12))
#
#     lum_veloc_name = 'lum-veloc_csmTrue_' + SN_name
#     lum_veloc_path = os.path.join(results_dir, lum_veloc_name)
#     lum_veloc_normalized_name = 'lum-veloc-normalized_csmTrue_' + SN_name
#     lum_veloc_normalized_path = os.path.join(results_dir, lum_veloc_normalized_name)
#
#     plot_snec_fits.plot_result_fit(lum_path, 'lum', axs[0, 0])
#     plot_snec_fits.plot_result_fit(lum_veloc_path, 'lum', axs[0, 1])
#     plot_snec_fits.plot_result_fit(lum_veloc_normalized_path, 'lum', axs[0, 2])
#     plot_snec_fits.plot_result_fit(lum_path, 'veloc', axs[1, 0])
#     plot_snec_fits.plot_result_fit(lum_veloc_path, 'veloc', axs[1, 1])
#     plot_snec_fits.plot_result_fit(lum_veloc_normalized_path, 'veloc', axs[1, 2])
#
#     axs[0, 0].set_ylabel('Log bolometric luminosity (erg/s)', fontsize=14)
#     axs[1, 0].set_ylabel('Expansion velocity (km/s)', fontsize=14)
#     axs[1, 0].set_xlabel('Rest-frame days from discovery', fontsize=14)
#     axs[1, 1].set_xlabel('Rest-frame days from discovery', fontsize=14)
#     axs[1, 2].set_xlabel('Rest-frame days from discovery', fontsize=14)
#
#     plt.tight_layout()
#
#     fig.savefig(os.path.join(output_dir, SN_name + '_lum_veloc_comparison.png'))
#
#     return fig



'''
for martinez in martinez_truth:
    for result_path in res_paths:
        res_name = result_path.split(sep='/')[-1]
        print(res_name)
        for n_step_str in n_steps_str:
            for sn_name in sn_names:
                if sn_name in result_path and n_step_str in result_path:
                    # check and if needed, make dirs in for the figure outputs
                    nstep_output_path = os.path.join(output_dir, n_step_str)
                    if not os.path.exists(nstep_output_path):
                        os.mkdir(nstep_output_path)
                    sn_name_output_path = os.path.join(output_dir, n_step_str, sn_name)
                    if not os.path.exists(sn_name_output_path):
                        os.mkdir(sn_name_output_path)
                    res_output_path = os.path.join(output_dir, n_step_str, sn_name, res_name)
                    if not os.path.exists(res_output_path):
                        os.mkdir(res_output_path)
                    # choose which plotting code to use
                    if 'csmFalse' in result_path:
                        if martinez:
                            res_martinez_output_path = os.path.join(output_dir, n_step_str, sn_name, res_name,
                                                                    'martinez_comparison')
                            if not os.path.exists(res_martinez_output_path):
                                os.mkdir(res_martinez_output_path)
                            plot_snec_fits_martinez_and_data_csmFalse.loopy_plot(result_path, res_martinez_output_path)
                        else:
                            plot_snec_fits_csmFalse.loopy_plot(result_path, res_output_path)
                    else:
                        if martinez:
                            res_martinez_output_path = os.path.join(output_dir, n_step_str, sn_name, res_name,
                                                                    'martinez_comparison')
                            if not os.path.exists(res_martinez_output_path):
                                os.mkdir(res_martinez_output_path)
                            plot_snec_fits_martinez_and_data.loopy_plot(result_path, res_martinez_output_path)
                        else:
                            plot_snec_fits.loopy_plot(result_path, res_output_path)

'''