import os
import sys
import numpy as np
import tensorflow as tf

import fastfiz as f
# from fastfiz import rules as r
game_type = 1 #'GT_EIGHTBALL'

def parse_gamestate(s_gamestate):
	s_splt = s_gamestate.split(' ')
	s_splt = [float(x) for x in s_splt]
	gameType = s_splt[0]
	turnType = s_splt[1]
	timeLeft = s_splt[2]
	timeLeft_opp = s_splt[3]
	curPlayer_started = s_splt[4]
	num_balls = s_splt[5]
	balls_arr = [s_splt[6+i*5:11+i*5] for i in range((balls_arr - 5) // 5)]
	return

def gen_shot_params(pseed):
	return


## Define shot params
## Simulate a shot
## Start a new frame
#Gamestate.rack()
# gamestate = f.GameState.RackedState(1)
# string_gamestate = gamestate.toString()


#GS = f.GameState_RackedState()
#shot = f.ShotParams(_a, _b, _theta, _phi, _v)
#shot.processShot(eventList, gameshot);


#if __name__=="__main__":
    
