from fastfiz import fastfiz

import cv2
import numpy as np
import os
import itertools

# Colors:
bg_color=(15,150,15)
color_steps=5
colors_array=np.linspace(0, 255, color_steps) 
# Remove green colored ball later
ball_colors=list(itertools.product(colors_array, colors_array, colors_array))
ball_colors.insert(0, (255, 255, 255))
ball_radius=2

def render_balls(dims, ball_positions, ball_params=None, output_filename='positions.png'):
	img_height, img_width, margin = dims
	img = np.zeros((img_width + 2*margin, img_height + 2*margin, 3), np.uint8)
	for x in range(margin, img_width + margin):
		for y in range(margin, img_height + margin):
			img[x, y, :] = bg_color[:]
	
	for b_id, b_pos in enumerate(ball_positions):
		cv2.circle(img, tuple(map(lambda x, y: x + y, b_pos, (margin, margin))), ball_radius, ball_colors[b_id], 6)

	cv2.imshow('img', img)
	cv2.waitKey(0)
	#img.save(output_filename, "PNG")

def get_ball_positions():
## TODO : Caller to fastfiz API

def get_game_config():
## TODO

if __name__=="__main__":
	img_width = 480
	img_height = 240
	margin = 20
	ball_positions = [(10, 10), (50, 50), (60, 60), (80, 60), (95, 10)]	

	render_balls([img_width, img_height, margin], ball_positions)
