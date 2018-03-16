#ifndef _FastFiz_h
#define _FastFiz_h

#ifdef SWIG
%include exception.i
%include "std_string.i"
%{
#include "FastFiz.h"
    %}
#endif /* SWIG */

#include <math.h>
#include <cassert>
#include <cstring>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <algorithm>
#include <set>
#include <string>
#include <list>
#include <vector>
#include <gsl/gsl_rng.h>

    //TODO: add a setprecision for double print out

    /**
     * The namespace in which all of the physics and rules components for FastFiz reside.  
     */
    namespace Pool {

        using namespace std;

        class BadShotException;
        class ShotParams;
        class Point;
        class Vector;
        class Ball;
        class TableState;
        class Table;
        class Event;
        class Shot; 
        class Utils;
#ifndef SWIG

        /**
         * The exception thrown when a Shot cannot be completed.
         * It can be thrown when TableState::executeShot() or TableState::getFirstBallHit() is called. 
         * The Type cannot be changed after construction.
         */
        class BadShotException : public exception {
            public:
                enum Type { OK,                 /**< The Shot was okay. */
                    MISSED_BALL,        /**< The Shot missed the CUE ball. */ 
                    INVALID_CUE_PARAMS, /**< The ShotParams for the Shot in the current TableState were out of range or not physically possible. */ 
                    POOLFIZ_ERROR,      /**< FastFiz encountered some sort of error. */
                    UNKNOWN_ERROR       /**< An unknown error. */ 
                };

                /**
                 * Constructs a BadShotException initialized with the given Type. 
                 * There is no default constructor.
                 */
                BadShotException (Type type) : _type(type) {}

                /**
                 * Returns the Type of the BadShotException. 
                 */
                Type getType() const {return _type;}

                /**
                 * Returns the Type of the BadShotException as a string.
                 */
                string getTypeString() const;

                /**
                 * Returns a const char* pointing to a c string about what type of Exception this is.
                 */
                virtual const char* what() const throw() {
                    return "Bad Shot Exception";
                }
            private:
                const Type _type;
        };
#endif /* SWIG */

        /**
         * A struct of the parameters for a CUE_STRIKE.  It is passed to TableState::executeShot() and can be retrieved from CueStrikeEvent. 
         * @see operator<<(ostream &out, const ShotParams &rhs) Writes human-readable output for ShotParams objects to the output stream.
         */
        struct ShotParams {
            double a;     /**< The x-coordinate of the cue stick (right/left english) on the CUE ball in mm. */
            double b;     /**< The y-coordinate of the cue stick (top/bottom english) on the CUE ball in mm. */
            double theta; /**< The elevation of the cue stick in degrees. */
            double phi;   /**< The azumith angle (heading) of the cue stick in degrees. */
            double v;     /**< The velocity of the cue stick before impact in m/s (max is 4.5 m/s). */

            /**
             * Constructs a ShotParams object with all the fields initialized to 0.0.
             */
            ShotParams() : a(0.0), b(0.0), theta(0.0), phi(0.0), v(0.0) {}

            /**
             * Constructs a ShotParams object with all the fields initialized to the given values.
             * @param _a     The x-coordinate of the cue stick (right/left english) on the CUE ball in mm. 
             * @param _b     The y-coordinate of the cue stick (top/bottom english) on the CUE ball in mm. 
             * @param _theta The elevation of the cue stick in degrees. 
             * @param _phi   The azumith angle (heading) of the cue stick in degrees. 
             * @param _v     The velocity of the cue stick before impact in m/s (max is 4.5 m/s). 
             */
            ShotParams(double _a, double _b, double _theta, double _phi, double _v) :
                a(_a), b(_b), theta(_theta), phi(_phi), v(_v) {}

            /**
             * Copy constructor.
             */
            ShotParams(const ShotParams &rhs) : 
                a(rhs.a), b(rhs.b), theta(rhs.theta), phi(rhs.phi), v(rhs.v) {}
        };

#ifndef SWIG
        /**
         * Writes human-readable output for the ShotParams object to the output stream.
         */
        ostream& operator<<(ostream &out, const ShotParams &rhs); 
#endif /* SWIG */

        /** 
         * A two-dimensional Point class. It has methods for a few simple point operations. 
         * @see operator<<(ostream &out, const Point &rhs) Writes human-readable output for Point objects to the output stream.
         */
        struct Point
        {
            double x; /**< The x-coordinate of the Point in meters. */
            double y; /**< The y-coordinate of the Point in meters. */

            /**
             * Constructs a Point object with x and y initialized to 0.0.
             */
            Point() : x(0.0), y(0.0) {}

            /**
             * Constructs a Point object with x and y initialized to the given values.
             */
            Point(double xx, double yy) : x(xx), y(yy) {}

            /**
             * Copy constructor.
             */
            Point(const Point &rhs) : x(rhs.x), y(rhs.y) {}

#ifndef SWIG
            /**
             * Calculates the rotation of a Point object about the origin.  The angle of rotation is counterclockwise.
             * @param cos_phi The cos of the angle to be rotated by. 
             * @param sin_phi The sin of the anlge to be rotated by.
             * @return A new Point object that is a copy of this Point object rotated about the origin by the angle 
             *         represented by the cos and sin values.
             */
            Point rotate(double cos_phi, double sin_phi) const;

            /**
             * Creates a Vector object representation of this Point object.
             * @return A new Vector object that is a copy of this Point object with the z field initialized to 0.0.
             */
            Vector to_v() const;

            /**
             * Calculates the sum of two Points. 
             * @param p2 The other Point object to be added. 
             * @return A new Point object that has the respective x and y components of the other Point object
             *         added to this Point object's components.
             */
            Point operator+(const Point& p2) const;

            /**
             * Returns a new Point object that has the respective x and y components of the other Point object
             * subtracted from this Point object's components.
             */
            Point operator-(const Point& p2) const;

            /**
             * Writes a machine-readable string representation of this Point object (a sequence of space-separated 
             * doubles followed by a trailing space) to the stream provided.  
             * It can later be interpreted by fromStream() of fromString().  
             */
            void toStream(ostream &out) const;

            /**
             * Reads a machine-readable string representation of a Point (a sequence of space-separated doubles)
             * written by toStream() or toString() from the stream provided and assigns the values to this Point object.
             */
            void fromStream(istream &in);

#endif /* ! SWIG */

            /**
             * Returns a string object that is a machine-readable string representation of this Point object
             * (a sequence of space-separated doubles followed by a trailing space).  
             * It can later be interpreted by fromString() or fromStream().
             */
            string toString() const;

            /**
             * Takes a machine-readable string representation of a Point (a sequence of space-separated doubles) 
             * created by toString() or toStream() and assigns the values to this Point object.  
             */
            void fromString(const string &s);
        };

#ifndef SWIG
        /**
         * Writes human-readable output for the Point object to the output stream.
         */
        ostream& operator<<(ostream &out, const Point &rhs);

#endif /* !SWIG */

        /** 
         * A three-dimensional Vector class. It has methods for various vector operations. 
         * @see operator<<(ostream &out, const Vector &rhs) Writes human-readable output for Vector objects to the output stream.
         */ 
        struct Vector
        {
            double x; /**< The x-coordinate of the Vector in meters. */
            double y; /**< The y-coordinate of the Vector in meters. */
            double z; /**< The z-coordinate of the Vector in meters. */

            /**
             * Constructs a Vector object with x, y and z initialized to 0.0.
             */
            Vector() : x(0.0), y(0.0), z(0.0) {}

            /**
             * Constructs a Point object with x, y and z initialized to the given values.
             */
            Vector(double xx, double yy, double zz) : x(xx), y(yy), z(zz) {}

            /**
             * Copy constructor.
             */
            Vector(const Vector &rhs) : x(rhs.x), y(rhs.y), z(rhs.z) {}

#ifndef SWIG
            /**
             * Returns a new Vector object that has the respective x, y, and z components of this Vector object
             * multiplied with the coefficient provided. 
             */
            Vector operator*(double d) const;

            /**
             * Multiplies the x,y,z coordinates of the vector by the coefficient provided.
             */
            Vector& operator*=(double d);

            /**
             * Returns a new Vector object that has the respective x, y, and z components of the other Vector object
             * subtracted from this Vector object's components.
             */
            Vector operator-(const Vector &fizV) const;

            /**
             * Returns a new Vector object that has the respective x, y, and z components of the other Vector object
             * added to this Vector object's components.
             */
            Vector operator+(const Vector &fizV) const;

            /**
             * Returns a new Vector object that is a copy of this Vector object rotated about the origin.
             * The parameter is in radians.  The angle is counterclockwise.
             */
            Vector rotateRad(double rad) const; //takes radians
            /**
             * Returns a new Vector object that is a copy of this Vector object rotated about the origin.
             * The parameter is in degrees.  The angle is counterclockwise.
             */
            Vector rotateDeg(double phi) const; //takes degrees
            /**
             * Returns a new Vector object that is a copy of this Vector object rotated about the origin.
             * The parameters are the cos and sin of the angle to be rotated by.  The angle is counterclockwise.
             */
            Vector rotate(double cos_phi, double sin_phi) const;
            /**
             * Returns dot-product of two vectors.
             */
            double dot(const Vector& v2) const;
            /**
             * Returns cross-product of two vectors.
             */
            Vector cross(const Vector& v2) const;
            /**
             * Returns the length (magnitude) of a vector.
             */
            double mag() const;
            /**
             * Returns the vector divided by its magnitude.
             */
            Vector norm() const;
            /**
             * Returns a point with same x,y as the vector.
             */
            Point to_p() const;

            void toStream(ostream &out) const;
            void fromStream(istream &in);
#endif /* ! SWIG */
            string toString() const;
            void fromString(const string &s);
        };
#ifndef SWIG
        ostream& operator<<(ostream &out, const Vector &rhs); //doesn't need to be friend
#endif /* ! SWIG */

#ifndef SWIG
        /** General utility functions */
        class Utils {
            public:
                static const double EPSILON = 1.0E-11;              /**< double comparison threshold */
                static const double VELOCITY_EPSILON = 1E-10;        /**< threshold velocity for balls to be considered stationary. */
                //Angle ranges can START_ZERO [0,2pi) or START_NEG_HALF_CIRCLE [-pi,pi)
                //This can apply to degrees or radians
                enum ANGLE_RANGE { START_ZERO, START_NEG_HALF_CIRCLE};
                enum ANGLE_UNIT {DEGREES,RADIANS};

                /** direction in the plane from p1 to p2 */
                static double angle(Point p1,Point p2,ANGLE_RANGE range,ANGLE_UNIT unit);
                /** converts radians to degrees */
                static double toDegrees(double angleRad);
                /** converts an angle in radians to be on the specified range, either [0,2*pi) or [-pi,pi) */
                static double normalizeRadianAngleRange(double angle, ANGLE_RANGE range);

                /** A random number generator instance */
                static gsl_rng* rng();

                /** Equality test for doubles. */
                static bool fequal(double a, double b);
                /** Less than test for doubles. */
                static bool fless(double a, double b);
                /** Greater than test for doubles. */
                static bool fgreater(double a, double b);
                /** Greater than or equal test for doubles. */
                static bool fgreaterequal(double a, double b);
                /** Less than or equal test for doubles. */
                static bool flessequal(double a, double b);
                /** Check if velocity absolute value should be considered zero. */
                static bool vzero(double a);    
            private:
                static gsl_rng* _rng;
        };
#endif /* ! SWIG */

        /** This class represents the state and position of a single pool ball.
         *  It is mostly used internally for execution of the physics, and as part of a TableState object.
         */
        class Ball
        {
            public:
#ifndef SWIG
                static const double MU_BALL_BALL = 0.01;            /**< coefficient of ball-ball friction */
                static const double MU_CUETIP_BALL = 0.7;           /**< coefficient of cue tip-ball friction */
                static const double BALL_COEFF_REST_POS = 2.0;      /**< coefficient of restitution for balls */
                static const double BALL_COEFF_REST_NEG = 0.0;      /**< coefficient of restitution for balls */
                static const double BALL_MASS = 163.01;             /**< 5.57 oz. ball mass [g] */
                static const double BALL_RADIUS = 0.028575;         /**< 2.25" ball radius [m] */
#endif /* SWIG */

                /** The present physical state of the ball. */
                enum State { NOTINPLAY, /**< not on the table */
                    STATIONARY, /**< in play but not moving */
                    SPINNING, /**< not moving translationally, but spinning in place about the vertical axis */
                    SLIDING, /**< sliding across the surface of the table, not necessarily in a straight trajectory */
                    ROLLING, /**< rolling across the surface of the table, necessarily in a straight trajectory */
                    POCKETED_SW, /**< pocketed in the SW pocket */
                    POCKETED_W,  /**< pocketed in the W pocket */
                    POCKETED_NW, /**< pocketed in the NW pocket */
                    POCKETED_NE, /**< pocketed in the NE pocket */
                    POCKETED_E,  /**< pocketed in the E pocket */
                    POCKETED_SE, /**< pocketed in the SE pocket */
                    SLIDING_SPINNING, /**< Transition between SLIDING and SPINNING */
                    ROLLING_SPINNING, /**< Transition between ROLLING and SPINNING */
                    UNKNOWN_STATE /**< in an unknown state */};

                /** The type (number) of the ball. */
                enum Type { CUE, ONE, TWO, THREE, FOUR, FIVE, SIX, SEVEN, EIGHT, 
                    NINE, TEN, ELEVEN, TWELVE, THIRTEEN, FOURTEEN,
                    FIFTEEN, UNKNOWN_ID };

                /** Create a new ball */
                Ball() : radius(BALL_RADIUS), r(), state(UNKNOWN_STATE), type(UNKNOWN_ID), v(), w() {}

                /** Create a new ball with specified number */
                Ball(Type _type) : radius(BALL_RADIUS), r(), state(UNKNOWN_STATE), type(_type), v(), w() {}
#ifndef SWIG
                /** Create a new ball with specified number and state */
                Ball(Type _type, State _state) : radius(BALL_RADIUS), r(), state(_state), type(_type), v(), w() {}
                /** Create a new ball with specified number, state, and position */
                Ball(Type _type, State _state, Point _r) : radius(BALL_RADIUS), r(_r), state(_state), type(_type), v(), w() {}
                /** Create a new ball with specified number, state, and position */
                Ball(Type _type, State _state, double x, double y) : radius(BALL_RADIUS), r(x,y), state(_state), type(_type), v(), w() {}
#endif /* ! SWIG */

                /** Copy constructor */
                Ball(const Ball &rhs) : radius(rhs.radius), r(rhs.r), state(rhs.state), type(rhs.type), v(rhs.v), w(rhs.w) {}
#ifndef SWIG
                /** Assignment operator */
                const Ball& operator= (const Ball &rhs);
#endif /* ! SWIG */

                /** Returns the radius of the ball in meters. */
                double getRadius() const {return radius;}
                /** Returns the ID (number) of the ball */
                Type getID() const {return type;}
                /** Returns the name of the ball (CUE, ONE, etc.) as a string */
                string getIDString() const; 
                /** Returns the physical state of the ball */
                State getState() const {return state;}
                /** Returns the physical state of the ball as a string */
                string getStateString() const;
                /** Returns the x,y position of the ball */
                const Point& getPos() const {return r;}
                /** Returns the x,y,z velocity of the ball */
                const Vector& getVelocity() const {return v;}
                /** Returns the x,y,z spin of the ball */
                const Vector& getSpin() const {return w;}

                /** Set the number of the ball */
                void setID(Type t) {type = t;}
                /** Set the position of the ball */
                void setPos(Point pos) {r = pos;}
                /** Set the velocity of the ball */
                void setVelocity(Vector vel) {v = vel;}
                /** Set the spin of the ball */
                void setSpin(Vector spin) {w = spin;}
                /** Set the physical state of the ball */
                void setState(State s) {state = s;}
                /** Returns true iff the ball is on the table. */
                bool isInPlay() const {return (state == STATIONARY || state == SPINNING || state == SLIDING || state == ROLLING);}
                /** Returns true iff the ball is in a pocket. */
                bool isPocketed() const {return (state == POCKETED_SW || state ==  POCKETED_W || 
                        state == POCKETED_NW || state == POCKETED_NE || 
                        state == POCKETED_E || state == POCKETED_SE);}
#ifndef SWIG
                /** Returns true if the ball occupies the same space as the other ball */
                bool overlaps(const Ball &other,double epsilon=Utils::EPSILON) const;

                void moveAway(Ball &other, double epsilon);

                /** Returns the distance between the centers of this ball and the other ball. 
                 * Note that a positive distance does not indicate no overlap.
                 * No overlap is ensured if the distance is greater than the sum of the radiuses of the balls.
                 */
                double dist(const Ball &other) const;
                /** Add small gaussian noise to ball's position based on dither parameter passed.
                 * This method is used for racking and spotting noise, and has nothing to do with shot execution noise.
                 */
                void addNoise(double dither);
#endif /* ! SWIG */
                /** Update balls physical state based on position, velocity, and spin.
                 *  Used internally for physics calculations.
                 */
                void updateState(bool VERBOSE = false);

#ifndef SWIG
                /** Print a machine-readable representation of the ball's type, radius and position.
                 * Only writes location information, not velocity and spin.
                 */
                void toStream(ostream &out) const;
                /** Read a machine-readable representation of the ball's type, radius and position.
                 * Should be used on data written by toStream()/toString().
                 */
                void fromStream(istream &in);
#endif /* ! SWIG */
                /** Returns a string including a machine-readable representation of the ball's type, radius and position.
                 * Only outputs location information, not velocity and spin.
                 */
                string toString() const;
                /** Read a machine-readable representation of the ball's type, radius and position from string.
                 * Should be used on data created by toString()/toStream().
                 */
                void fromString(const string &s);
#ifndef SWIG

            private:
                double radius;
                Point r;
                State state;
                Type type;
                Vector v;
                Vector w;
#endif /* ! SWIG */
        };

#ifndef SWIG
        ostream& operator<<(ostream &out, const Ball &rhs); //doesn't need to be a friend
#endif /* !SWIG */

        /** Class representing physical properties of a pool table. This does not include the position of the balls. 
         * In general, this object doesn't need to be instansiated more than once.
         * This class also provides some static utility functions and constants related to the pool table.
         */
        class Table
        {
            public:
                /* Physics constants */
                static const double g = 9.81;                       /**< Gravitational constant.*/
                static const double MU_SLIDING = 0.2;		/**< coefficient of sliding friction */
                static const double MU_ROLLING = 0.015;		/**< coefficient of rolling friction */
                static const double MU_SPINNING = 0.044;		/**< coefficient of spinning friction */

                /* Default table parameters */
                static const double TABLE_LENGTH = 2.236;		/**< table length [m] */
                static const double TABLE_WIDTH = 1.116;		/**< table width [m] */
                static const double CORNER_POCKET_WIDTH = 0.11;	/**< corner pocket width [m] */
                static const double SIDE_POCKET_WIDTH = 0.12;	/**< side pocket width [m] */
                static const double RAIL_HEIGHT = 0.040005;		/**< rail height [m] */
                static const double CUE_LENGTH = 1.45;		/**< cue length [m] */
                static const double RAIL_VEL_DAMPING_X = 0.6;	/**< damping factor for velocity component parallel to rail */
                static const double RAIL_VEL_DAMPING_Y = 0.9;	/**< damping factor for velocity component perpendicular to rail */
                static const double RAIL_SPIN_DAMPING = 0.1;	/**< damping factor for angular velocity component */
                static const double RAIL_VEL_ANGLE_ADJ = 0.0;	/**< angle adjustment factor for velocity vector */
                static const double RAIL_ZSPIN_ANGLE_ADJ = 0.0;	/**< angle adjustment factor for vertical component of angular velocity vector */

                //used in CueStrikeEvent::doHandle
                static const double CUE_MASS = 600.0;
                static const double I = 0.4*160.0*0.028575*0.028575*0.995473; /**< moment of inertia; to match poolfiz */

                /** A rail or a pocket. */
                enum BoundaryId { SW_POCKET, SW_RAIL, W_POCKET, NW_RAIL, NW_POCKET, N_RAIL, NE_POCKET,
                    NE_RAIL, E_POCKET, SE_RAIL, SE_POCKET, S_RAIL, UNKNOWN_BOUNDARY };

                enum Pocket { SW, W, NW, NE, E, SE, UNKNOWN_POCKET };

                /** Constructs a table with the default parameters */
                Table() : 
                    _muS(MU_SLIDING), 
                    _muR(MU_ROLLING), 
                    _muSp(MU_SPINNING), 
                    _railHeight(RAIL_HEIGHT), 
                    _cueLength(CUE_LENGTH), 
                    _railVelDampingX(RAIL_VEL_DAMPING_X), 
                    _railVelDampingY(RAIL_VEL_DAMPING_Y), 
                    _railSpinDamping(RAIL_SPIN_DAMPING), 
                    _railZSpinAngleAdj(RAIL_ZSPIN_ANGLE_ADJ), 
                    _railVelAngleAdj(RAIL_VEL_ANGLE_ADJ) 
            {
                calculateTableDimensions( TABLE_LENGTH, TABLE_WIDTH, CORNER_POCKET_WIDTH, SIDE_POCKET_WIDTH);
            }

                /** Constructs a table with the given parameters and automatically sets the pockets.
                 *	@param length  the length of the table in metres
                 *	@param width the width of the table in metres
                 *	@param cornerPocketWidth the width of corner pockets (horn to horn) in metres 
                 *	@param sidePocketWidth the width of the side pockets (horn to horn) in metres
                 *	@param muS the coefficient of sliding friction (dimensionless)
                 *	@param muR the coefficient of rolling friction (dimensionless)
                 *	@param muSp the coefficient of spinning friction (dimensionless)
                 *	@param railHeight the height of the top of the rail above the table in metres
                 *	@param cueLength the length of the cue in metres 
                 *	@param railVelDampingX velocity damping factor of the banks (X)
                 *  @param railVelDampingY velocity damping factor of the banks (Y)
                 *	@param railSpinDamping spin damping factor of the banks
                 *	@param railZSpinAngleAdj z-spin angle of deflection factor of the banks
                 *	@param railVelAngleAdj velocity deflection factor of the banks */
                Table( double length, 
                        double width, 
                        double cornerPocketWidth, 
                        double sidePocketWidth, 
                        double muS = MU_SLIDING,
                        double muR = MU_ROLLING,
                        double muSp = MU_SPINNING,
                        double railHeight = RAIL_HEIGHT,
                        double cueLength = CUE_LENGTH,
                        double railVelDampingX = RAIL_VEL_DAMPING_X,
                        double railVelDampingY = RAIL_VEL_DAMPING_Y,
                        double railSpinDamping = RAIL_SPIN_DAMPING,
                        double railZSpinAngleAdj = RAIL_ZSPIN_ANGLE_ADJ,
                        double railVelAngleAdj = RAIL_VEL_ANGLE_ADJ);

                /** Copy constructor. */
                Table( const Table &rhs );

#ifndef SWIG
                /** Assignment operator */
                const Table& operator= ( const Table &rhs );
#endif /* SWIG */

                /* Convenience methods for the table parameters */

                /** Returns the length of the table (in metres). */
                double getLength() const {return _length;}

                /** Returns the width of the table (in metres). */
                double getWidth() const {return _width;}

                /** Returns the y-coordinate of the headstring of the table (in metres) */
                double getHeadString() const {return _headString;}

                /** Returns the x-y coordinate of the footspot of the table (in metres) */
                const Point& getFootSpot() const {return _footSpot;}

                /** Sets the length of the cue stick. 
                 *	@param length the length of the cue in metres. */
                void setCueLength(double length) {_cueLength = length;}

                /** Returns the length of the cue stick. */
                double getCueLength() const {return _cueLength;}

                /** Sets the height of the rails around the table.
                 *	@param height the height of the rails in metres. */
                void setRailHeight(double height) {_railHeight = height;}

                /** Returns the height of the rails around the table. */
                double getRailHeight() const {return _railHeight;}

                /** Sets the coefficient of sliding friction.
                 *	@param mu the coefficient of sliding friction (between 0 and 1). */
                void setMuSliding(double mu) {_muS = mu;} 

                /** Returns the coefficient of sliding friction. */
                double getMuSliding() const {return _muS;}

                /** Sets the coefficient of rolling friction.
                 *	@param mu the coefficient of rolling friction (between 0 and 1). */
                void setMuRolling(double mu) {_muR = mu;}

                /** Returns the coefficient of rolling friction. */
                double getMuRolling() const {return _muR;}

                /** Sets the coefficient of spinning friction.
                 *	@param mu the coefficient of spinning friction (between 0 and 1). */
                void setMuSpinning(double mu) {_muSp = mu;}

                /** Returns the coefficient of spinning friction. */
                double getMuSpinning() const {return _muSp;}

                /** Returns the aiming point of the given pocket (in metres). */
                const Point & getPocketCenter(Pocket pocket) const;

                /** Returns the right corner coordinates of the given pocket (in metres), looking 
                 *	at the pocket from the center of the table. */
                const Point & getPocketRight(Pocket pocket) const;

                /** Returns the left corner coordinates of the given pocket (in metres), looking 
                 *	at the pocket from the center of the table. */
                const Point & getPocketLeft(Pocket pocket) const;

                /** Contains a single default table object whose reference is returned when called */
                static const Table& defaultTable() {
                    static const Table t = Table();
                    return t;
                }

                /** Returns the pocketed ball state corresponding to the given pocket (by Pocket) */
                static Ball::State stateFromPocket(Pocket pocket);
                /** Returns the same pocket in Pocket form rather than BoundaryId form */
                static Pocket pocketFromBndId(BoundaryId bnd);
                /** Returns the same pocket in BoundaryId form rather than Pocket form */
                static BoundaryId bndIdFromPocket(Pocket pocket);
                /** returns a string of the boundary name */
                static string boundaryName(BoundaryId boundary);
                /** returns a string of the pocket name */
                static string pocketName(Pocket pocket);
#ifndef SWIG

            private:
                double _muS;
                double _muR;
                double _muSp;
                double _railHeight;
                double _cueLength;
                double _railVelDampingX;
                double _railVelDampingY;
                double _railSpinDamping;
                double _railZSpinAngleAdj;
                double _railVelAngleAdj;

                double _length;
                double _width;
                double _headString;
                Point _footSpot;
                Point _SWpocketLeft;
                Point _SWpocketRight;
                Point _SWpocketCenter;
                Point _WpocketLeft;
                Point _WpocketRight;
                Point _WpocketCenter;
                Point _NWpocketLeft;
                Point _NWpocketRight;
                Point _NWpocketCenter;
                Point _NEpocketLeft;
                Point _NEpocketRight;
                Point _NEpocketCenter;
                Point _EpocketLeft;
                Point _EpocketRight;
                Point _EpocketCenter;
                Point _SEpocketLeft;
                Point _SEpocketRight;
                Point _SEpocketCenter;

                void calculateTableDimensions( double length, double width, double cornerPocketWidth, double sidePocketWidth );
#endif /* ! SWIG */
        };

        /** Base class of the Event heirarchy. Indicates a notable event in the execution of a shot.
         * 
         * Subclasses of this class indicate specific events.
         */
        class Event
        {
#ifndef SWIG
            friend ostream& operator<<(ostream &out, const Event &rhs);
#endif /* ! SWIG */
            public:
            /** Indicates the type of event.*/
            enum Type { NO_EVENT, /**< Nothing happened (only used internally) */
                STATE_CHANGE, /**< A ball transitioned to a different motion state */
                BALL_COLLISION, /**< two balls collided */
                RAIL_COLLISION, /**< a ball collided with a rail */
                POCKETED, /**< a ball was pocketed */
                CUE_STRIKE, /**< the cue ball was struck by the cue stick successfully */
                MISCUE,  /**< the cue ball was struck by the cue stick unsuccessfully */
                UNKNOWN_EVENT  /**< an unknown event */
            };

            /** Create a new base event.
             * \param time Simulation time event took place.
             * \param b id of the ball affected by the event.
             */
            Event(double time, Ball::Type b) : _time(time), _ball1(b), _ball1Data(NULL) {}

            /** Returns the event's simulation time index */
            double getTime() const {return _time;}

            /** Returns the id of first (possibly only) ball involved in the event. */
            Ball::Type getBall1() const {return _ball1;}
            /** Returns the physical information (location, velocity,spin) of the first (possibly only) ball involved in the event.
             * May cause a null pointer exception if the information is not available (if event was not yet handled).
             */
            Ball& getBall1Data() {return *_ball1Data;}

            /** Sorts Events by time */
            bool operator<(const Event &other) const;
            /** Event comparison function, calls operator< */
            static bool eventCmp(const Event* event1, const Event* event2);

            /** Returns a human-readble string representation of the event. */
            string toString() const;

            /** Returns the type of the event */
            virtual Type getType() const {return UNKNOWN_EVENT;}
            /** Returns the type of the event as a string */
            virtual string getTypeString() const {return "Unknown Event";}
            /** Returns the id of the second ball involved in the event, if applicable */
            virtual Ball::Type getBall2() const {return Ball::UNKNOWN_ID;}
            /** Returns the physical information (location, velocity,spin) of the second ball involved in the event.
             * If no second ball is involved, will return the first ball information.
             * May cause a null pointer exception if the information is not available (if event was not yet handled).
             */
            virtual Ball& getBall2Data() {return *_ball1Data;}
            /** Returns true if this event and the other event affect a shared ball. */
            virtual bool relatedTo(const Event &other) const;
            /** Returns true if this event affects specified ball. */
            virtual bool involvesBall(Ball::Type b) const;
            /** Destructor. */
            virtual ~Event() {delete _ball1Data;}

            /** Used within physics code to apply event to a table state, modifying ball states, positions, and properties 
             * Also updates ball information within event.
             */
            void handle(TableState &ts, bool VERBOSE = false) {doHandle(ts, VERBOSE); copyBalls(ts);}
            protected:
            double _time;
            Ball::Type _ball1;
            Ball* _ball1Data;

            virtual void doHandle(TableState &ts, bool VERBOSE) const = 0;
            virtual void copyBalls(TableState &ts);
#ifndef SWIG
            virtual ostream& dump(ostream &out) const;
            private:
            //copy and assignment
            Event(const Event &rhs); //: _time(rhs._time), _ball1(rhs._ball1) {}
        const Event& operator=(const Event & rhs); //{if (this != &rhs) {_time=rhs._time; _ball1=rhs._ball1;}; return *this;}
#endif /* ! SWIG */
        };

        /** An event invloving a physical transition of state, e.g between SLIDING and ROLLING.
         * In general these events are important for graphical simulation of the shot, and for trajectory analysis.
         */
        class StateChangeEvent : public Event {
            public:
                StateChangeEvent(double time, Ball::Type b) : Event(time, b) {}

                virtual Type getType() const {return STATE_CHANGE;}
                virtual string getTypeString() const {return "State Change";}
            protected:
                virtual void doHandle(TableState &ts, bool VERBOSE) const;
#ifndef SWIG
            private:
                //copy and assignment
                StateChangeEvent(const StateChangeEvent &rhs); // : Event(rhs) {}
        StateChangeEvent& operator=(const StateChangeEvent & rhs); //{if (this != &rhs) Event::operator=(rhs); return *this;} 
#endif /* ! SWIG */
        };

        /** An event involving two balls colliding. 
         * This is the only event where Ball2 is defined.
         */
        class BallCollisionEvent : public Event {
            public:
                BallCollisionEvent (double time, Ball::Type b1, Ball::Type b2) : Event(time,b1), _ball2(b2), _ball2Data(NULL) {}

                virtual Type getType() const {return BALL_COLLISION;}
                virtual string getTypeString() const {return "Ball Collision";}
                virtual bool relatedTo(const Event &other) const;
                virtual bool involvesBall(Ball::Type b) const;
                virtual ~BallCollisionEvent() {delete _ball2Data;}

                virtual Ball::Type getBall2() const {return _ball2;}
                virtual Ball& getBall2Data() {return *_ball2Data;}
            protected:
                Ball::Type _ball2;
                Ball* _ball2Data;

                virtual void copyBalls(TableState &ts);
                virtual void doHandle(TableState &ts, bool VERBOSE) const;
#ifndef SWIG
                virtual ostream& dump(ostream &out) const;
            private:
                //copy and assignment
                BallCollisionEvent (const BallCollisionEvent &rhs); //: Event(rhs), _ball2(rhs._ball2) {}
        BallCollisionEvent & operator=(const BallCollisionEvent& rhs);
#endif /* ! SWIG */
        };

        /** An event involving a ball colliding with a rail. */
        class RailCollisionEvent : public Event {
            public:
                RailCollisionEvent (double time, Ball::Type b, Table::BoundaryId rail) 
                    : Event(time,b), _rail(rail) {}

                virtual Type getType() const {return RAIL_COLLISION;}
                virtual string getTypeString() const {return "Rail Collision";}

                /** Returns the ID of the rail collided with. */
                Table::BoundaryId getRail() const {return _rail;}
            protected:
                Table::BoundaryId _rail;

                virtual void doHandle(TableState &ts, bool VERBOSE) const;
#ifndef SWIG
                virtual ostream& dump(ostream &out) const;
            private:
                //copy and assignment
                RailCollisionEvent (const RailCollisionEvent &rhs); //: Event(rhs), _rail(rhs._rail) {}
        RailCollisionEvent & operator=(const RailCollisionEvent & rhs);
#endif /* ! SWIG */
        };

        /** An event involving a ball entering a pocket. */
        class PocketedEvent : public Event {
            public:
                PocketedEvent (double time, Ball::Type b, Table::Pocket pocket) 
                    : Event(time,b), _pocket(pocket) {}

                virtual Type getType() const {return POCKETED;}
                virtual string getTypeString() const {return "Ball Pocketed";}

                /** Returns the ID of the pocket. */
                Table::Pocket getPocket() const {return _pocket;} 
            protected:
                Table::Pocket _pocket;

                virtual void doHandle(TableState &ts, bool VERBOSE) const;
#ifndef SWIG
                virtual ostream& dump(ostream &out) const;
            private:
                //copy and assignment
                PocketedEvent (const PocketedEvent &rhs); //: Event(rhs), _pocket(rhs._pocket) {}
        PocketedEvent& operator=(const PocketedEvent & rhs);
#endif /* ! SWIG */
        };

        /** The initial event in any shot -- cue ball struck by the cue stick. */
        class CueStrikeEvent : public Event {
            public:
                //Default constructor
                /** Constructor with shot parameters, assumes zero time */
                CueStrikeEvent (const ShotParams &params) 
                    : Event(0.0, Ball::CUE), _params(params) {}
                /** Constructor with shot parameters, and time */
                CueStrikeEvent (double time, const ShotParams &params) 
                    : Event(time, Ball::CUE), _params(params) {}
                /** Constructor with shot parameters, ball, and time */
                CueStrikeEvent (double time, Ball::Type &b, const ShotParams &params) 
                    : Event(time,b), _params(params) {}

                virtual Type getType() const {return CUE_STRIKE;}
                virtual string getTypeString() const {return "Cue Strike";}

                /** Return shot parameters */
                const ShotParams& getParams() const {return _params;} 
            protected:
#ifndef SWIG
                ShotParams _params;
#endif /* SWIG */
                virtual void doHandle(TableState &ts, bool VERBOSE) const;
#ifndef SWIG
                virtual ostream& dump(ostream &out) const;
            private:
                //copy and assignment
                CueStrikeEvent (const CueStrikeEvent &rhs); //: Event(rhs), _params(rhs._params) {}
        CueStrikeEvent& operator=(const CueStrikeEvent& rhs);
#endif /* ! SWIG */
        };


        /** Detailed result of a simulation of a shot.
         * Public objects of this class can only be created by TableState::executeShot()
         * Constructing a Shot object starts the physics simulation.
         * The private methods in this class implement the main part of the physics simulation.
         */
        class Shot
        {
            public:    
                friend class TableState;

                /** Gets the list of events generated during the shot.
                 * The vector is sorted by time, and all relevant ball information is available in the events.
                 */
                const vector<Event*>& getEventList() const {return shotEvents;}
                /** Returns the amount of time (in seconds) balls will take to settle after executing the shot. */
                double getDuration() const;

                /** Destroy a shot object. 
                 * Deleting the Shot object is caller's responsibility.
                 */
                ~Shot();
            private:
                //Uncopiable
                Shot(const Shot &rhs);
                Shot& operator=(const Shot &rhs);

#ifndef SWIG
                vector<Event*> shotEvents;

                //Can throw BadShotException from simulateShot
                Shot(TableState &state,const ShotParams &sp, bool fullSim /* false = get first ball hit only */,
                        bool verbose = false, bool errors=false) : shotEvents() {
                    simulateShot(state,sp,fullSim,verbose,errors);
                }

                //Can throw BadShotException
                void simulateShot(TableState &state,
                        const ShotParams &sp,
                        bool fullSim,
                        bool verbose,
                        bool errors);

                static bool VERBOSE;
                static bool CHECK_ERRORS;
                static void updateTime(TableState &state, double oldTime, double newTime);
                static void updateBall(const Table &table, Ball &ball, double oldTime, double newTime);
                static double updateSpinning(double w_z, double t, double mu_sp, double R, bool isSliding);

                static void removeRelatedEvents(list<Event*> &events, Event *lastEvent);
                static void addRelatedFutureEvents(TableState &state, double curTime, list<Event*> &futureEvents, Ball::Type b1, Ball::Type b2);
                static Event* nextTransitionEvent(const Table &table, Ball &ball, double curTime);
                static Event* nextCollisionEvent(const Table &table, Ball &ball1, Ball &ball2, double curTime);
                static void addBoundaryEvents(const Table &table, Ball &ball, double curTime, list<Event*> &futureEvents);
                static double calcEventTime(int numRoots, double root1, double root2, double curTime);

                static int solveQuartic(double roots[], double a0, double a1, double a2, double a3, double a4);
                static const double NEAR_FUTURE_EPSILON = 1E-8;
                static double leastPositiveRealRoot(double roots[],double epsilon=NEAR_FUTURE_EPSILON);


                /*
                   static void handleEvent(TableState &state, Event *event);
                   static void handleStateChange(Ball &ball);
                   static void handleBallCollision(Ball &ball1, Ball &ball2); 
                   static void handleRailCollision(Ball &ball, Table::BoundaryId rail);
                   static void handlePocketed(Ball &ball, Table::Pocket pocket);
                   static void handleCueStrike(Ball &ball, const ShotParams &sp);
                   */

                //static void addFutureEvents(TableState &state, double curTime, list<Event*> &futureEvents);

                //static const double BALL_MASS = 160.0; //From Marlow, in line with BCA rules
                //static const double g = GSL_CONST_MKSA_GRAV_ACCEL;
                //static const double R = 0.028575; //Ball radius //Now taken from Ball::BALL_RADIUS
                //static const int NUM_BALLS = 16; //Now taken from TableState::getNumBalls()
#endif /* ! SWIG */
        };
#ifndef SWIG
        /** Print a human-readable representation of the shot to a stream */
        ostream& operator<<(ostream &out, Shot &rhs); //doesn't need to be friend
#endif /* ! SWIG */

        /** The physical state of balls on a table.
         *  Includes position and velocity information but no rules related information.
         */
        class TableState
        {
            public:
                /* Thresholds */
                static const double MAX_VELOCITY = 10;		/**< maximum velocity allowed [m/s]  */
                static const double MIN_THETA = 0.0;		/**< minimum theta allowed [degrees] */
                static const double MAX_THETA = 70.0;		/**< maximum theta allowed [degrees] */

                /** This integer is used within the bitmask integer returned by Table::isPhysicallyPossible()
                 *	and Table::isValidBallPlacement() (as well as the fizG versions of those classes). It 
                 *	indicates that all shot parameters are valid; the actual integer returned by these methods 
                 *	is simply the OR'ed combinations of the all the possibilities from the list under the heading
                 *	'Variables' in the HTML documentation or listed below in the header file; you can AND the result 
                 *	returned by these methods with any of the possibilities from the list below to test if that 
                 *	particular possibility is true.  */
                static const int OK_PRECONDITION = 0;

                /** This integer is used within the bitmask integer returned by Table::isPhysicallyPossible()
                 *	and Table::isValidBallPlacement() (as well as the fizG versions of those classes).  It 
                 *	indicates that shot parameter 'a' is invalid or not physically possible;
                 *	the actual integer returned by these methods 
                 *	is simply the OR'ed combinations of the all the possibilities from the list under the heading
                 *	'Variables' in the HTML documentation or listed below in the header file; you can AND the result 
                 *	returned by these methods with any of the possibilities from the list below to test if that 
                 *	particular possibility is true.  */
                static const int BAD_A_VAL = 1; 	

                /** This integer is used within the bitmask integer returned by Table::isPhysicallyPossible()
                 *	and Table::isValidBallPlacement() (as well as the fizG versions of those classes).  It 
                 *	indicates that shot parameter 'b' is invalid or not physically possible;
                 *	the actual integer returned by these methods 
                 *	is simply the OR'ed combinations of the all the possibilities from the list under the heading
                 *	'Variables' in the HTML documentation or listed below in the header file; you can AND the result 
                 *	returned by these methods with any of the possibilities from the list below to test if that 
                 *	particular possibility is true.  */
                static const int BAD_B_VAL = 2;

                /** This integer is used within the bitmask integer returned by Table::isPhysicallyPossible()
                 *	and Table::isValidBallPlacement() (as well as the fizG versions of those classes).  It 
                 *	indicates that shot parameter 'theta' is invalid or not physically possible;
                 *	the actual integer returned by these methods 
                 *	is simply the OR'ed combinations of the all the possibilities from the list under the heading
                 *	'Variables' in the HTML documentation or listed below in the header file; you can AND the result 
                 *	returned by these methods with any of the possibilities from the list below to test if that 
                 *	particular possibility is true.  */
                static const int BAD_THETA_VAL = 4;

                /** This integer is used within the bitmask integer returned by Table::isPhysicallyPossible()
                 *	and Table::isValidBallPlacement() (as well as the fizG versions of those classes).  It 
                 *	indicates that shot parameter 'phi' is invalid or not physically possible;
                 *	the actual integer returned by these methods 
                 *	is simply the OR'ed combinations of the all the possibilities from the list under the heading
                 *	'Variables' in the HTML documentation or listed below in the header file; you can AND the result 
                 *	returned by these methods with any of the possibilities from the list below to test if that 
                 *	particular possibility is true.  */
                static const int BAD_PHI_VAL = 8;

                /** This integer is used within the bitmask integer returned by Table::isPhysicallyPossible()
                 *	and Table::isValidBallPlacement() (as well as the fizG versions of those classes).  It 
                 *	indicates that shot parameter 'V' is invalid or not physically possible;
                 *	the actual integer returned by these methods 
                 *	is simply the OR'ed combinations of the all the possibilities from the list under the heading
                 *	'Variables' in the HTML documentation or listed below in the header file; you can AND the result 
                 *	returned by these methods with any of the possibilities from the list below to test if that 
                 *	particular possibility is true.  */
                static const int BAD_V_VAL = 16; 	

                /** This integer is used within the bitmask integer returned by Table::isPhysicallyPossible()
                 *	and Table::isValidBallPlacement() (as well as the fizG versions of those classes).  It 
                 *	indicates that the cue ball placement x-coordinate is out of range;
                 *	the actual integer returned by these methods 
                 *	is simply the OR'ed combinations of the all the possibilities from the list under the heading
                 *	'Variables' in the HTML documentation or listed below in the header file; you can AND the result 
                 *	returned by these methods with any of the possibilities from the list below to test if that 
                 *	particular possibility is true.  */
                static const int BAD_X_VAL = 32;	

                /** This integer is used within the bitmask integer returned by Table::isPhysicallyPossible()
                 *	and Table::isValidBallPlacement() (as well as the fizG versions of those classes).  It 
                 *	indicates that the cue ball placement y-coordinate is out of range;
                 *	the actual integer returned by these methods 
                 *	is simply the OR'ed combinations of the all the possibilities from the list under the heading
                 *	'Variables' in the HTML documentation or listed below in the header file; you can AND the result 
                 *	returned by these methods with any of the possibilities from the list below to test if that 
                 *	particular possibility is true.  */
                static const int BAD_Y_VAL = 64;

                /** This integer is used within the bitmask integer returned by Table::isPhysicallyPossible()
                 *	and Table::isValidBallPlacement() (as well as the fizG versions of those classes).  It 
                 *	indicates that the cue will collide with a ball before striking the cue ball;
                 *	the actual integer returned by these methods 
                 *	is simply the OR'ed combinations of the all the possibilities from the list under the heading
                 *	'Variables' in the HTML documentation or listed below in the header file; you can AND the result 
                 *	returned by these methods with any of the possibilities from the list below to test if that 
                 *	particular possibility is true.  */
                static const int CUE_STICK_COLLISION = 128;

                /** This integer is used within the bitmask integer returned by Table::isPhysicallyPossible()
                 *	and Table::isValidBallPlacement() (as well as the fizG versions of those classes).  It 
                 *	indicates that the ball placement x-y coordinates result in an overlap with another ball on table;
                 *	the actual integer returned by these methods 
                 *	is simply the OR'ed combinations of the all the possibilities from the list under the heading
                 *	'Variables' in the HTML documentation or listed below in the header file; you can AND the result 
                 *	returned by these methods with any of the possibilities from the list below to test if that 
                 *	particular possibility is true.  */
                static const int BALL_OVERLAP = 256;

                /** Create an empty table state for given table (or default table if not specified) */
                TableState(const Table &table=Table::defaultTable()) : _table(table), balls() {}
                /** Copy constructor */
                TableState(const TableState &rhs) : _table(rhs._table), balls(rhs.balls) {}
#ifndef SWIG
                /** Assignment operator */
                TableState& operator=(const TableState &rhs);
#endif /* ! SWIG */

                /** Return the number of different balls used in representing this state */
                int getNumBalls() const {return balls.size();}
#ifndef SWIG
                /** Iterator for start of ball vector */
                vector<Ball>::iterator getBegin() {return balls.begin();}
                /** Iterator for end of ball vector */
                vector<Ball>::iterator getEnd() {return balls.end();}
                /** Const Iterator for start of ball vector */
                vector<Ball>::const_iterator getBegin() const {return balls.begin();}
                /** Const Iterator for end of ball vector */
                vector<Ball>::const_iterator getEnd() const {return balls.end();}

                void fixOverlap(const bool VERBOSE);
#endif /* ! SWIG */
                /** Modify specified ball id with specified ball information, adding if needed */
                void setBall(Ball &b);
                /** Set a ball's state and position, adding if needed */
                void setBall(Ball::Type btype, Ball::State state, Point r);
                /** Set a ball's state and position, adding if needed */
                void setBall(Ball::Type btype, Ball::State state, double x, double y);
                /** Spot a pocketed ball at or near the foot spot, with dithering if needed */
                void spotBall(Ball::Type btype, double dither = DITHER);
                /** Get ball information based on id. If ball does not exist it will be added with unknown information. 
                 * Returned ball may be modified and will affect the table state.
                 */
                Ball& getBall(Ball::Type btype);
#ifndef SWIG
                /** Get ball information based on id. If ball does not exist it will be added with unknown information. */
                const Ball& getBall(Ball::Type btype) const;
#endif /* SWIG */

                /** Return the underlying table object */
                const Table& getTable() const {return _table;};

                /* Validity checkers */

                /** Returns BAD_X_VAL or BAD_Y_VAL if any ball that is in play is off the table or
                 *	BALL_OVERLAP if any ball is overlapping with another ball.  
                 *	Used to check cue ball placement validity in the case of Ball In Hand, for example.
                 *	@param VERBOSE Print debugging information */
                int isValidBallPlacement(bool VERBOSE = false) const;

                /** Checks physical validity of the given shot parameters and table state.  Used before calling
                 *	executeShot() to double-check that an error condition will not result due to invalid
                 *	parameters. 
                 *  @param shotParams Shot parameters
                 *	@param VERBOSE Print debugging information
                 */
                int isPhysicallyPossible(const ShotParams &shotParams, bool VERBOSE=false) const;

                /** Add racking noise to table state.
                 * This has nothing to do with shot execution noise.
                 */
                void addNoise(double dither);
#ifndef SWIG
                /** Randomly position all balls on the table. */
                void randomize();
#endif /* ! SWIG */

#ifdef SWIG
                %exception {
                    try {
                        $function
                    } catch (Pool::BadShotException& ex) {
                        string msg="Bad Shot: "+ex.getTypeString();
                        SWIG_exception(SWIG_ValueError,msg.c_str());
                    }
                }

                %newobject executeShot;
#endif /* SWIG */
                /** Execute given shot on the table state. 
                 * The TableState object will be modified to reflect the final state of the shot.
                 * If you need to retain the initial state, copy the object before calling.
                 * Noise is not added to the shot prior to execution.
                 * @param sp The shot parameters to execute.
                 * @param verbose Print a lot of debugging information to cerr.
                 * @param errors Detect internal errors (slower execution)
                 * @exception BadShotException will be thrown in case of bad parameters or errors in execution.
                 */
                Shot* executeShot(const ShotParams &sp, bool verbose=false, bool errors=false) {
                    return new Shot(*this, sp, true, verbose, errors);
                }

#ifdef SWIG
                %exception;
#endif /* SWIG */

                //Returns the first ball that will be hit by the cue ball by the given shot.
                //If the cue ball doesn't hit another ball, returns CUE.
                //Can throw BadShotException from creating Shot
                /** Return the first ball hit by the cue ball while executing shot.
                 * This function will simulate the shot up until the first ball is hit by 
                 * the cue ball and then return the ball type, or UNKNOWN_ID if no ball is hit.
                 */
                Ball::Type getFirstBallHit(const ShotParams &sp);

#ifndef SWIG
                /** Writes a machine-readable representation of the table state to a stream. */
                void toStream(ostream &out) const;
                /** Reads a machine-readable representation of the table state from a stream. 
                 * Should be used on output of toStream() or toString().
                 */
                void fromStream(istream &in);
#endif /* ! SWIG */
                /** Returns a machine-readable representation of the table state as a string*/
                string toString() const;
                /** Reads a machine-readable representation of the table state from a string. 
                 * Should be used on output of toStream() or toString().
                 */
                void fromString(const string &s);
#ifndef SWIG

            private:
                const Table &_table;
                vector<Ball> balls;

                static const double DITHER = 0.00005;                 
                static const double EPSILON_B = 0.001;		/**< small value added to the b shot parameter for threshold calculations [mm] */
                static const double EPSILON_THETA = 0.001;		/**< small value added to the theta shot parameter for threshold calculations [degrees] */
                static const int MAX_BALLS_ON_TABLE = 100;          /**< Maximum balls allowed on the table at any time. */

                static int numLineSphereIntersections(Vector &p1, Vector &p2, Vector &p3, double rad, double& root1, double& root2);
                int findOverlap(Ball::Type ball) const;  
#endif /* ! SWIG */
        };
#ifndef SWIG
        ostream& operator<<(ostream &out, TableState &rhs); //doesn't need to be friend
#endif /* ! SWIG */

        /** Return a string identifying the version and build information of the library */
        string getFastFizVersion();
#ifndef SWIG
        /** Return a string identifying the version and build information of the library.
         *  @deprecated Use getFastFizVersion().
         */
        inline string getPoolfizVersion() {return getFastFizVersion();}
#endif /* ! SWIG */

#ifndef SWIG
        /** Rack table state for Eight-ball. Debugging function.
         * @deprecated Use the GameState framework to generate a racked state.
         */
        void rack(TableState& ts);
        /** Dump event list to standard error. Debugging function.
        */
        void printEvents(list<Event*> events);
        /** Print overlapping balls in table state. Debugging function.
        */
        void printOverlap(TableState& ts);
#else /* SWIG */
        %newobject getTestState;
#endif /* SWIG */
        /** Return a racked state for testing purposes. Debugging function. */
        TableState* getTestState();
#ifdef SWIG
        %newobject getTestShotParams();
#endif /* SWIG */
        /** Return sample shot parameters that work with the test state. */
        ShotParams* getTestShotParams();

    }//namespace Pool

#ifdef SWIG

%include "std_vector.i"

namespace std{
    %template(EventVector) vector<Pool::Event*>;
}
#endif /* SWIG */

#endif //_FastFiz_h
