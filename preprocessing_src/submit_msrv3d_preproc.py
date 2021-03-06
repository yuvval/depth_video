#!/usr/bin/python
import submit_jobs_matlab
import scipy.io
from numpy import *
import itertools

    

prepr_params = dict();

dir_names = submit_jobs_matlab.runBash('ls /storage/yuvval/data/video/MSRV3D/').split()

prepr_params['eigen.clgtv'] = """struct('depth_method', 'eigen', 'opflow_method', 'CLGTV', 'sample_interval', {iv})"""

frame_smp_interval = 1

run_prms = "preprocess_wrapper('msrv3d', {vid_name}, {prepr_params})"

dir_names = [dir for dir in dir_names if 'Indoor' in dir]
num_vids_in_dirs = []
vid_names = []
for dir in dir_names:
    info = scipy.io.loadmat('/storage/yuvval/data/video/MSRV3D/' + dir + '/info.mat')
    num_vids_in_dir = shape(info['clips'])[0]
    vid_names += ["struct('dir', '%s', 'idx', %d)"%(dir, k+1) for k in range(num_vids_in_dir)]
 

Jobs = []
# Jobs += [("ctx06", run_prms, 'approaching_toward.avi', prepr_params['fayao.clgtv'].format(iv=frame_smp_interval))]

for k in range(len(vid_names)):
    Jobs += [("baffin", run_prms, vid_names[k], prepr_params['eigen.clgtv'].format(iv=frame_smp_interval))]
# for k in range(15, 30):
#     Jobs += [("ctx07", run_prms, video_names[k], prepr_params['fayao.clgtv'].format(iv=frame_smp_interval))]
# for k in range(30, 45):
#     Jobs += [("ctx08", run_prms, video_names[k], prepr_params['fayao.clgtv'].format(iv=frame_smp_interval))]
# for k in range(45, 60):
#     Jobs += [("ctx10", run_prms, video_names[k], prepr_params['fayao.clgtv'].format(iv=frame_smp_interval))]
# for k in range(60, 75):
#     Jobs += [("ctx11", run_prms, video_names[k], prepr_params['fayao.clgtv'].format(iv=frame_smp_interval))]
#
# for k in range(75, 95):
#     Jobs += [("ctx17", run_prms, video_names[k], prepr_params['fayao.clgtv'].format(iv=frame_smp_interval))]


kill_flag = False
#jobs_to_kill = [job[0] for job in Jobs] # killing all jobs in Jobs
#jobs_to_kill = ['ctx%02d'%ctx_num for ctx_num in range(4,19)] # killing jobs according to cortex machine numbers


if kill_flag:
    submit_jobs_matlab.kill_jobs(jobs_to_kill)
else:
    submit_jobs_matlab.submit_jobs(Jobs, ignore_git=True)



# bag1
# basketball1
# basketball2
# basketball2.2
# basketballnew
# bdog_occ2
# bear_back
# bear_change
# bird1.1_no
# bird3.1_no
# book_move1
# book_turn
# book_turn2
# box_no_occ
# br_occ_0
# br_occ1
# br_occ_turn0
# cafe_occ1
# cc_occ1
# cf_difficult
# cf_no_occ
# cf_occ2
# cf_occ3
# child_no2
# computerbar1
# computerBar2
# cup_book
# dog_no_1
# dog_occ_2
# dog_occ_3
# express1_occ
# express2_occ
# express3_static_occ
# face_move1
# face_occ2
# face_occ3
# face_turn2
# flower_red_occ
# gre_book
# hand_no_occ
# hand_occ
# library2.1_occ
# library2.2_occ
# mouse_no1
# new_ex_no_occ
# new_ex_occ1
# new_ex_occ2
# new_ex_occ3
# new_ex_occ5_long
# new_ex_occ6
# new_ex_occ7.1
# new_student_center1
# new_student_center2
# new_student_center3
# new_student_center4
# new_student_center_no_occ
# new_ye_no_occ
# new_ye_occ
# one_book_move
# rose1.2
# static_sign1
# studentcenter2.1
# studentcenter3.1
# studentcenter3.2
# three_people
# toy_car_no
# toy_car_occ
# toy_green_occ
# toy_mo_occ
# toy_no
# toy_no_occ
# toy_wg_no_occ
# toy_wg_occ
# toy_wg_occ1
# toy_yellow_no
# tracking4
# tracking7.1
# two_book
# two_dog_occ1
# two_people_1.1
# two_people_1.2
# two_people_1.3
# walking_no_occ
# walking_occ1
# walking_occ_long
# wdog_no1
# wdog_occ3
# wr_no
# wr_no1
# wr_occ2
# wuguiTwo_no
# zball_no1
# zball_no2
# zball_no3
# zballpat_no1
