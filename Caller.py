import os
import sys
import numpy as np
# import tensorflow as tf

sys.path.append('./fastfiz/fastfiz/')
# sys.path.append('./fastfiz/')

from fastfiz import fz as f
# from fastfiz import rules as r

game_type = 1 #'GT_EIGHTBALL'


"""
Parse the Gamestate from string to ball positions.
GameState format: 
			GameType, TurnType, timeleft, timeleft_opp, curplayerstarted, numballs
			foreach ball : radius, state, type, Point.x, Point.y 
			1 1 (unknown)
"""
def parse_gamestate(s_gamestate):
	s_splt = s_gamestate.split(' ')
	s_splt = s_splt[:-3]
	s_splt = [float(x) for x in s_splt]
	gameType = s_splt[0]
	turnType = s_splt[1]
	timeLeft = s_splt[2]
	timeLeft_opp = s_splt[3]
	curPlayer_started = s_splt[4]
	num_balls = s_splt[5]
	balls_arr = [s_splt[6+i*5:11+i*5] for i in range((len(s_splt) - 5) // 5)]
	return balls_arr

def gen_shot_params(pseed):
	return


if __name__=="__main__":
	## Define shot params
	## Simulate a shot
	## Start a new frame
	#Gamestate.rack()
	gamestate = f.GameState.RackedState(1)
	string_gamestate = gamestate.toString()
	ball_arr = parse_gamestate(string_gamestate)
	print(ball_arr)
	# shot = f.ShotParams(5., 5., 10., 30., 2.)
	# Gshot = f.GameShot()
	# Gshot.ball = 1
	# Gshot.params = shot
	# Gshot.pocket = 0
	# gamestate.executeShot(Gshot)
	# string_gamestate_after = gamestate.toString()
	# ball_arr_after = parse_gamestate(string_gamestate_after)
	# print(np.array(ball_arr_after)-np.array(ball_arr))
	# shot.processShot(eventList, gameshot);
