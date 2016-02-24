""" Some utility methods such as:
Printing and executing command lines.
Copying a directory tree.
"""
# Author: Bishesh Khanal <bishesh.khanal@inria.fr>
#
#

import subprocess
import shutil

def print_and_execute(command, print_cmd=True):
    """prints the command by default and executes it in the shell"""
    if print_cmd:
        print command + '\n'
    subprocess.call(command, shell=True)

def update_or_execute_cmd(cmd, cmd_n, execute):
    '''
    Either update or execute cmd.
    If execute is True executes cmd_n command, returns cmd.
    otherwise returns cmd + cmd_n.
    '''
    if execute:
        print_and_execute(cmd_n)
        return cmd
    else:
        return cmd + cmd_n
def get_lines_as_list(fname, delim):
    '''
    Return a list that contains all the lines present in the file fname
    except those line that starts with delim.
    '''
    with open(fname, 'r') as fil:
        return [line.strip() for line in fil if not line.startswith(delim)]

def copy_dir_tree(src, dest):
    """ Copy entire directory tree rooted at the src directory.
    Uses shutil interface."""
    try:
        shutil.copytree(src, dest)
    # Directories are the same
    except shutil.Error as ex:
        print('Directory not copied. Error: %s' % ex)
    # Any error saying that the directory doesn't exist
    except OSError as ex:
        print('Directory not copied. Error: %s' % ex)


def rem_num_of_lines(in_filename, start_string):
    """
    Returns the number of lines that is remaining in the given file starting
    from the first appearance of the given start string.
    ---------
    PARAMETERS:
    <in_filename> Input filename

    <start_string> Start string to be search in <in_file>.

    ---------
    RETURNS:
    <lines_rem> Number of lines remaining in <in_file> after the first
    appearance of <start_string>, includes the line where <start_string> first
    appeared in the count.
    """
    lines_rem = 0
    start_string_found = False
    with open(in_filename) as in_file:
        for line in in_file:
            if not start_string_found:
                if start_string in line:
                    start_string_found = True
            if start_string_found is True:
                lines_rem += 1
    return lines_rem


from popen2 import popen2
def qsub_job(name, queue, walltime, procs, mem, dest_dir, cmd, add_atrb=None):
    """
    Adapted from:
    https://waterprogramming.wordpress.com/2013/01/24/pbs-job-submission-with-python/
    Submit jobs with qsub.
    To add if email address option desired:
    #PBS -M your_email@address
    #PBS -m abe  # (a = abort, b = begin, e = end)
    ---------
    PARAMETERS
    Parameters up to <mem> are self-explanatory and they set corresponding pbs
    script options.

    <dest_dir> Directory where job output and error files are to be written.

    <add_atrb> sting that will be used with -W option of qsub. These are
    additional_atributes for e.g. useful to make the job dependent on other jobs
    .

    <cmd> All the desired commands you want to launch in the job.
    -----------
    RETURNS
    <job_id> The job id provided by the system.
    """
    # get qsub options and actual sets of commands to run from the job
    job_string = """#!/bin/bash
#PBS -N %s
#PBS -q %s
#PBS -l walltime=%s
#PBS -l %s
#PBS -l %s
#PBS -d %s
""" % (name, queue, walltime, procs, mem, dest_dir)
    if add_atrb is None:
        job_string = job_string + ('cd $PBS_O_WORKDIR \n %s\n' % cmd)
    else:
        job_string = job_string + ('#PBS -W %s \n cd $PBS_O_WORKDIR \n %s\n' %
                                   add_atrb, cmd)
    #Open a pipe to the qsub command:
    output, inpt = popen2('qsub')

    #Send job_string to qsub
    inpt.write(job_string)
    inpt.close()

    # Print the submitted job and the system response to the screen
    print(job_string)
    job_id = output.read()
    print(job_id)
    return job_id


def sophia_nef_pbs_setting():
    """ Some commands to custom set-up the inria sophia nef cluster. These are
    suggested set ups before your commands in a job script.
    """
    return """

# bind a mpi process to a cpu
export OMPI_MCA_mpi_paffinity_alone=1
export LD_LIBRARY_PATH=/opt/openmpi/current/lib64

"""

def oarsub_job(l, p, d, n, job, run_from_file=False):
    '''
    Submit an oarsub job. Single letter args corresponds to oarsub arg.
    -l l (resource), -p p (resource property), -d d (directory to launch
    command from.), -n n (job name).
    E.g. l corresponds to arguments to be given to -l of oarsub.
    If run_from_file, job is a script file name so uses -S, else inline string.
    '''
    cmd = 'oarsub'
    if l:
        cmd += (' -l %s' % (l))
    if p:
        cmd += (' -p %s' % (p))
    if d:
        cmd += (' -d %s' % (d))
    if n:
        cmd += (' -n %s' % (n))
    #cmd = 'oarsub -l %s -p %s -d %s -n %s ' % (l, p, d, n)
    #cmd = 'oarsub -l %s -d %s -n %s ' % (l, d, n)
    # job_prefx = ('export OMPI_MCA_mpi_paffinity_alone=1\n'
    #              'export LD_LIBRARY_PATH=/opt/openmpi/current/lib64\n')
    job_prefx = ''
    if run_from_file:
        cmd += ' -S %s' % (job)
    else:
        cmd += ' "%s%s"' % (job_prefx, job)
    #print cmd
    print_and_execute(cmd)
