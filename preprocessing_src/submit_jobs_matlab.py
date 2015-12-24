#!/usr/local/bin/python
import subprocess
import os
import time

def runBash(cmd):
    """ This function takes Bash commands and returns them """
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    out = p.stdout.read().strip()
    return out  #This is the stdout from the shell command
def seqlist(l):
	return ' '.join([str(x) for x in l])

def kill_jobs(Jobs):
    for job in Jobs: 
        cmd  = 'ssh ' + job + ' " killall MATLAB" '
        print cmd
        runBash(cmd)

def submit_jobs(Jobs, ignore_git  = False):
    """ running parallel matlab jobs on different machines 
    Jobs are indicated by the variable 'Jobs' . 
    Each record holds 3 fields: machine_address_string, matlab_function_to_call, a list of parameter to call with
    e.g.
    run_prms['calt+easy10.comet'] = "similarity_experiment('Caltech256+split=easy10cat', 'COMET+v2+hinge_logdet_Nfrob_PSD', 5, {hyp_prm_sweep}, {experiment_stage})"
    Jobs = ( 
            ("ctx01", run_prms["calt+easy10.comet"], "struct('step_size',[0.1 1 10 100 1000], 'logdet_fact',[0.01 0.1 0.5 1 5])", "'clean_junk_mutex'"),
            ("ctx01", run_prms["calt+easy10.comet"], "struct('step_size',[0.1 1 10 100 1000], 'logdet_fact',[0.01 0.1 0.5 1 5])", "'search_hyperparams'"),
            ("ctx04", run_prms["calt+easy10.comet"], "struct('step_size',[0.1 1 10 100 1000], 'logdet_fact',[0.01 0.1 0.5 1 5])", "'search_hyperparams'"),
        )
    """ 

    NumGitModified=int(runBash('git status | grep modified | wc -l'))
    NumGitNewUncommited=int(runBash("git status | grep  \"new file\" | wc -l"))
    NumGitUntrackedMatlab=int(runBash("git status | grep '.m$' | grep -v \"new file\" | wc -l"))

    if (NumGitUntrackedMatlab > 0):
        print "WARNING: Below MATLAB files are (git) untracked"
        print runBash("""git status | grep '.m$' | grep -v "new file" | grep -v "modified" """)


    if (NumGitModified > 0 or NumGitNewUncommited > 0):
        print "ERROR: Please (git) commit or revert the below modified files!!"
        print runBash("""git status | grep modified""")
        print runBash("""git status | grep "new file" """)
        if ignore_git:
            print 'Press any key to continue, OR CTRL+C to abort !!'
            txt = raw_input()
        else:
            exit()

    if (NumGitUntrackedMatlab > 0):
        print "Warning: some matlab files are untracked. type y to continue"
        txt = raw_input()
        if (txt != 'y'):
            print "JOBS WEREN'T EXECUTED!!!"
            exit()

    for k, job in enumerate(Jobs):
        matlab_script_txt = "addpath('../');\ninitdirs;\n{job[1]}".format(job=job).format(vid_name=job[2], prepr_params=job[3])
        matlab_script_fname = "%s_%d.m"%(job[0],k+1000)
        with open(matlab_script_fname, "w") as f:
            f.write(matlab_script_txt + "\n")

        cmd  = 'ssh ' + job[0] + ' " cd ~/depth_video/preprocessing_src/; source ../theano-env/bin/activate; ./run_bg_matlab_jobs_single_machine.sh ' + matlab_script_fname + '"'

        print cmd
        print runBash(cmd)
        if 'clean_junk_mutex' in matlab_script_txt:
            time.sleep(15)




