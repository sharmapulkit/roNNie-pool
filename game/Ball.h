#ifndef _BALL
#define _BALL

class Ball{
private:
	vect position;
	vect velocity;
	vect omega;
	float mass;
public:
	Ball();
	Ball(vect position, vect velocity, vect omega, float mass);

	//getters
	vect get_position();
	vect get_velocity();
	vect get_omega();
	
	//setters
	void set_position();
	void set_velocity();
	void set_omega();
};

#endif _BALL