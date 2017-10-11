#ifndef _VECT
#define _VECT

class vect{
private:
	float x;
	float y;
	float z;
public:
	//constructor
	vect();
	vect(float x, float y, float z);

	// getters
	float getX();
	float getY();
	float getZ();

	// setters
	void setX(float a);
	void setY(float a);
	void setZ(float a);

	// Utility
	float mag();
	void _add(vect a);	//adds vector a to this
	vect add(vect a);	//adds vector a and this, returns resultant
	void _sub(vect a);	//subtracts vector a from this
	vect sub(vect a);	//subtracts vector a from this returns resultant
	void _dot(vect a);	//dot product vector a to this 
	vect dot(vect a);	//dot product vector a to this returns resultant
	void _cross(vect a);	//cross product vector a to this returns resultant
	vect cross(vect a);		//cross product vector a to this returns resultant
	void _scalar(float a);	// multiplies this with scalar 
	vect scalar(float a); 	// multiplies this and scalar and returns resultant

	//Destructor
	~vect();
};

#endif