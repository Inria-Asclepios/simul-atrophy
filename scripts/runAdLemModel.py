#!/usr/bin/env python
import os
import os.path as op
import sys
import argparse as ag
import bish_utils as bu

def get_input_options():
    """ command line interface, get input options and interact with the user
    """
    #parse the command line inputs
    parser = ag.ArgumentParser()
    parser.add_argument(
        'patient', help='patient dir at $ADLEM_DIR/results/patients/ which '
        'contains all the inputs required by this script. Results will also be '
        ' written here.')
    parser.add_argument(
        'res_prefix', help='prefix to be added to all the output filenames')
    parser.add_argument('lame_paras', help='input lame parameters: '
                        'Format: "muTissue muCsf lambdaTissue lambdaCsf" '
                        'If both --mu_image and --in_dti used, lame_paras wont'
                        ' be used and those images will be used instead. ')
    parser.add_argument('atrophy', help='input atrophy file.')
    parser.add_argument('in_seg', help='input brain seg ')
    parser.add_argument(
        'in_img', help='valid input image that will be warped by the obtained'
        ' displacement fields from the model.')
    parser.add_argument('time_steps', help='(integer) '
                        'number of time steps you want to solve the system.')
    parser.add_argument('-m', '--mu_file', help='input mu image file: ')
    parser.add_argument('--use_dti', dest='use_dti', action='store_true',
                        help='If the option is provided, input DTI must be '
                        'given with the option -d')
    parser.add_argument('-d', '--in_dti', help='valid tensor image file present'
                        ' in patientID directory.')
    parser.add_argument('--wrt_press', action='store_true',
                        help='Write pressure map on the disk.')
    parser.add_argument('--wrt_force', action='store_true',
                        help='Write force on the the disk.')
    parser.add_argument('--wrt_residual', action='store_true',
                        help='Write residual on the disk.')
    parser.add_argument('--in_cluster', action='store_true',
                        help='If provided, launches as a job on the cluster.')
    parser.add_argument('--no_petsc_summary', action='store_true',
                        help='If provided, does not detail the petsc outputs.')

    ops = parser.parse_args()
    if ops.use_dti:
        if ops.in_dti is None:
            ops.in_dti = raw_input('Enter a DTI file (since you are using '
                                   '--use_dti.):  ')
    return ops


def get_petsc_options(no_detailed_output):
    """ Return petsc options. Change the code as required for different
    options you want.
    """
    ops = ('-pc_type fieldsplit -pc_fieldsplit_type schur '
           '-pc_fieldsplit_schur_precondition self -pc_fieldsplit_dm_splits 0 '
           '-pc_fieldsplit_0_fields 0,1,2 -pc_fieldsplit_1_fields 3 '
           '-fieldsplit_0_pc_type hypre')
    #monitor options
    # ops = ops + ('-fieldsplit_0_ksp_converged_reason -fieldsplit_0_ksp_max_it'
    #            ' 100 -ksp_converged_reason')
    if not no_detailed_output:
        ops = ops+(' -fieldsplit_1_ksp_converged_reason -ksp_converged_reason')
        # ops = ops + ('-fieldsplit_1_ksp_max_it 3 -ksp_max_it 3 -ksp_rtol 1.0e-8 '
        #              '-ksp_monitor_true_residual')
        ops = ops + (' -fieldsplit_1_ksp_monitor_true_residual '
                     '-ksp_monitor_true_residual -log_summary -ksp_view')
    return ops


def main():
    """
    Run the model adlem for a given patient, segmentation mask, atrophy map and
    [optional DTI]."""

    work_dir = os.getenv('ADLEM_DIR')
    if work_dir is None:
        raise ValueError('environment variable ADLEM_DIR not set.')
    target = op.join(work_dir, 'build/src/AdLemMain')
    ops = get_input_options()
    res_dir = op.join(work_dir, 'results/patients', ops.patient)
    in_img = op.join(res_dir, ops.in_img)
    atrophy = op.join(res_dir, ops.atrophy)
    in_seg = op.join(res_dir, ops.in_seg)
    res_path = res_dir+'/' #due to stupid limitation in my AdLem code
    # where it concatenates simply the -resPath argument with the file prefix.
    if ops.use_dti is True:
        use_dti = 'true'
        in_dti = op.join(res_dir, ops.in_dti)
    else:
        use_dti = 'false'
        in_dti = "dummyFileName"
    if ops.wrt_press is True:
        wrt_press = 'true'
    else:
        wrt_press = 'false'
    if ops.wrt_force is True:
        wrt_force = 'true'
    else:
        wrt_force = 'false'
    if ops.wrt_residual is True:
        wrt_residual = 'true'
    else:
        wrt_residual = 'false'
    petsc_ops = get_petsc_options(ops.no_petsc_summary)
    if ops.mu_file:
        mu_file = op.join(res_dir, ops.mu_file)
        mu_ops = '-parameters "%s" -muFile %s' % (ops.lame_paras, mu_file)
    else:
        mu_ops = '-parameters "%s"' % (ops.lame_paras)

    cmd = ('%s %s -atrophyFile %s -maskFile %s -imageFile %s -useTensorLambda '
           '%s -lambdaFile %s -numOfTimeSteps %s -resPath %s '
           '-resultsFilenamesPrefix %s -writePressure %s -writeForce %s '
           '-writeResidual %s %s'
           % (target, mu_ops, atrophy, in_seg, in_img, use_dti, in_dti,
              ops.time_steps, res_path, ops.res_prefix, wrt_press, wrt_force,
              wrt_residual, petsc_ops)) # check spaces in petsc_ops.
    # cmd = target + ' -atrophyFile ' + atrophy + ' -useTensorLambda ' + use_dti+\
    #       ' -lambdaFile ' + in_dti + ' -maskFile ' + in_seg + ' -imageFile ' + \
    #       in_img + ' -numOfTimeSteps ' + ops.time_steps + ' -resPath '+res_path\
    #       + ' -resultsFilenamesPrefix ' + ops.res_prefix + ' -writePressure ' +\
    #       wrt_press + ' -writeForce ' + wrt_force + ' -writeResidual ' + \
    #       wrt_residual + petsc_ops
    if ops.in_cluster is True:
        cluster_mpi = '/opt/openmpi-gcc/current/bin/mpiexec '
        cmd = bu.sophia_nef_pbs_setting() + cluster_mpi + cmd
        job_name = '%sSteps%s' % (ops.res_prefix, ops.time_steps)
        #print cmd
        bu.qsub_job(name=job_name, queue='parlong', walltime='12:00:00',
                    procs='nodes=3:xeon:ppn=20', mem='mem=180gb',
                    dest_dir=res_dir, cmd=cmd)
    else:
        #print cmd
        bu.print_and_execute(cmd)


if __name__ == "__main__":
    main()
