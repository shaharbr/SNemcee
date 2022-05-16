import os
from matplotlib import pyplot as plt
import plot_snec_fits
# import plot_snec_fits_martinez_and_data
import sys, getopt
import datetime


time_now = datetime.datetime.now().strftime("%Y-%m-%d_%H-%M-%S")

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


def plot_single(fig_type, model_path, ax):
    plot_snec_fits.plot_result_fit(model_path, fig_type, ax)
    ax.set_xlabel('Rest-frame days from discovery', fontsize=14)
    if fig_type == 'lum':
        ax.set_ylabel('Expansion velocity (km/s)', fontsize=14)
    if fig_type == 'veloc':
        ax.set_ylabel('Log bolometric luminosity (erg/s)', fontsize=14)
    if fig_type == 'mag':
        ax.set_ylabel('Absolute Magnitude', fontsize=14)
    plt.tight_layout()
    return ax


def composite_plot(SN_name, fig_type, fitting_type, csm, normalization, LumThreshold, results_dir, output_dir):
    fig_types = fig_type.split("-")
    model_name = SN_name + '_' + fitting_type + '_csm-' + csm + '_normalized' + str(normalization)
    if 'lum' in fitting_type:
        model_name += '_TreshLum' + str(LumThreshold)
    if 'mag' in fitting_type:
        model_name += '_TreshMagTrue'
    if 'veloc' in fitting_type:
        model_name += '_TreshVelocTrue'
    model_path = os.path.join(results_dir, model_name)
    num_subplots = len(fig_types)
    fig, axs = plt.subplots(num_subplots, figsize=(10, num_subplots*7))
    if num_subplots > 1:
        for i, f in enumerate(fig_types):
            plot_single(f, model_path, axs[i])
    else:
        plot_single(fig_type, model_path, axs)
    fig.savefig(os.path.join(output_dir, model_name + '_'+str(fig_type)+'_plot.png'))


def corner_plot(SN_name, fitting_type, csm, normalization, LumThreshold, results_dir, output_dir):
    model_name = SN_name + '_'+str(fitting_type)+'_csm-'+csm+'_normalized'+str(normalization)+'_TreshLum'+str(LumThreshold)
    model_path = os.path.join(results_dir, model_name)
    plot_snec_fits.overlay_corner_plot([model_path], output_dir,
                                       [model_name], model_name)


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
    fig, axs = plt.subplots(3, 2, figsize=(20, 12))

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


def lum_veloc_onestep_vs_twostep(SN_name, results_dir, output_dir):
    fig, axs = plt.subplots(2, 3, sharey='row', figsize=(20, 12))

    onestep_name = SN_name + '_lum-veloc_csm-with_normalizedFalse_TreshLumFalse_TreshVelocTrue'
    onestep_path = os.path.join(results_dir, onestep_name)
    twostep_priorNone_name = SN_name + '_lum_csm-twostep_normalizedFalse_TreshLumFalse'
    twostep_priorNone_path = os.path.join(results_dir, twostep_priorNone_name)
    twostep_priorTrue_name = SN_name + '_lum_csm-twostep-carryover_normalizedFalse_TreshLumFalse'
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
    SN_name = ''
    n_steps = 500
    type_fig = False
    output_dir = 'output_' + time_now
    fitting_type = False
    csm = False
    normalization = False
    LumThreshold = False

    arg_help = '{0} -S <SN name> -s <number of steps> -f <figure_type> -o <output directory> -t <fitting type> -c <csm> -n <normalization> -l <luminosity threshold>\n'\
               '\nargs:\n' \
               '-S SN name [required. for example: SN2017eaw]\n' \
               '-s number of steps = <int> [default: 500]\n' \
               '-f figure type = single fit-type plots: lum/veloc/mag [or any combination of those, separated by -] OR corner OR comparison plots: csm_comparison/lum-mag-veloc_comparison/lum-mag-veloc-Tthresh_comparison/lum-veloc-normalized/lum-veloc-twostep [required]\n' \
               '-o output directory name = <str> [default: output_<current time>]\n' \
               '-t fitting_type = lum, veloc, mag, lum-veloc, mag-veloc, lum-mag, lum-veloc-mag, combined [default: False] \n' \
               '-c csm = with, without, twostep, twostep-carryover [default: False] \n' \
               '-n normalization = True, False [default: False] \n' \
               '-l luminosity_threshold = True, False [default: False] \n' \
               ''.format(argv[0])
    try:
        opts, args = getopt.getopt(argv[1:], "hS:s:f:o:t:c:l:n:", ["help", "SN=", "steps=", "type of figure=", "output_dir=",
                                                                   "fitting type=", "csm=", "normalization=", "luminosity threshold="])
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
        elif opt in ("-f", "--figure_type"):
            type_fig = arg
        elif opt in ("-o", "--output_dir"):
            output_dir = arg
        elif opt in ("-t", "--fitting_type"):
            fitting_type = arg
        elif opt in ("-c", "--csm"):
            csm = arg
        elif opt in ("-n", "--normalization"):
            normalization = arg
        elif opt in ("-l", "--luminosity_threshold"):
            LumThreshold = arg
    print('SN name:', SN_name)
    print('steps:', n_steps)
    print('output directory:', output_dir)
    print('type of figure:', type_fig)
    print('fitting type:', fitting_type)
    print('csm:', csm)
    print('normalization:', normalization)
    print('luminosity threshold:', LumThreshold)

    res_dir = os.path.join('mcmc_results', output_dir, str(n_steps) + 'step')
    step_dir = os.path.join('figures', output_dir, str(n_steps) + 'step')
    if not os.path.exists(step_dir):
        if not os.path.exists(os.path.join('figures', output_dir)):
            os.mkdir(os.path.join('figures', output_dir))
        os.mkdir(step_dir)

    if 'comparison' not in type_fig:
        if type_fig == 'corner':
            corner_plot(SN_name, fitting_type, csm, normalization, LumThreshold, res_dir, step_dir)
        else:
            composite_plot(SN_name, type_fig, fitting_type, csm, normalization, LumThreshold, res_dir, step_dir)
    elif type_fig == 'csm_comparison':
        lum_wCSM_vs_woCSM(SN_name, res_dir, step_dir)
    elif type_fig == 'lum-mag-veloc_comparison':
        lum_veloc_vs_mag_veloc(SN_name,res_dir, step_dir)
    elif type_fig == 'lum-mag-veloc-Tthresh_comparison':
        lum_veloc_vs_mag_veloc(SN_name,res_dir, step_dir, LumTthresh=True)
    elif type_fig == 'lum-veloc-normalized_comparison':
        lum_vs_lum_veloc_vs_lum_veloc_normalized(SN_name,res_dir, step_dir)
    elif type_fig == 'lum-veloc-twostep_comparison':
        lum_veloc_onestep_vs_twostep(SN_name,res_dir, step_dir)


if __name__ == "__main__":
    main(sys.argv)


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