// 
#include "FastFiz.h"
#include "Stopwatch.h"

#include <gsl/gsl_sf_trig.h>
#include <gsl/gsl_poly.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_randist.h>

#define POOLFIZ_ERRORS
//#define ALL_EVENTS

#ifdef POOLFIZ_ERRORS
#define POOLFIZ_ERROR() {BadShotException _ex(BadShotException::POOLFIZ_ERROR); throw _ex;};
#else
#define POOLFIZ_ERROR() {}
#endif


//#include <gsl/gsl_complex.h>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_const_mksa.h>
//#include <gsl/gsl_complex_math.h>

/***********************************************************************************
 * Equations used in this physics simulator were taken from "An Event-based Pool   *
 * Physics Simulator" by Will Leckie and Michael Greenspan (the makers of Poolfiz) *
 * and from The Physics of Pocket Billiards by Wayland C. Marlow.                  *
 ***********************************************************************************/

namespace Pool {

    bool Shot::VERBOSE = false;
    bool Shot::CHECK_ERRORS = false;

    /*
       FILE *FastFiz::LOGFILE = NULL;

#ifndef NDEBUG
#include <signal.h>
#include "IOUtils.h"
static fizTableState dbgState;
static ShotParams dbgParams;

static void error(int signum) {
cerr << "ERROR CONDITION" << endl;
IOUtils::writeState(cerr,dbgState);
cerr << " " << dbgParams << endl;
}

#endif
*/
inline double Utils::toDegrees(double angleRad)
{
    return angleRad * 180 / M_PI;
}

inline bool Utils::vzero(double a) {
    return (fabs(a) < VELOCITY_EPSILON );
}

inline bool Utils::fequal(double a, double b) {
    return (fabs(a-b) < EPSILON );
}

inline bool Utils::fless(double a, double b) {
    return  (a < b && !fequal(a,b));
}

inline bool Utils::fgreater(double a, double b) {
    return (a > b && !fequal(a,b));
}

inline bool Utils::fgreaterequal(double a, double b) {
    return (a > b || fequal(a,b));
}

inline bool Utils::flessequal(double a, double b) {
    return (a < b || fequal(a,b));
}

//direction in the plane from p1 to p2
inline double Utils::angle(Point p1,Point p2,ANGLE_RANGE range,ANGLE_UNIT unit)
{
    double raw = atan2(p2.y-p1.y,p2.x-p1.x);
    double normed = normalizeRadianAngleRange(raw,range);
    if(unit == RADIANS)
    {
        return normed;
    }
    return toDegrees(normed);
}

double Utils::normalizeRadianAngleRange(double angle,ANGLE_RANGE range)
{
    //first, convert to [0,2*pi)
    double twoPi = 2 * M_PI;
    if(angle < 0)
    {
        double frac = angle / twoPi;
        double residual = frac - ((int)frac);
        angle = (1 + residual) * twoPi;
    }
    if(angle >= twoPi)
    {
        double frac = angle / twoPi;
        double residual = frac - ((int)frac);
        angle =  residual * twoPi;
    }

    if(START_NEG_HALF_CIRCLE == range)
    {
        if(angle >= M_PI)
        {
            angle -= twoPi;
        }
    }

    return angle;
}


string BadShotException::getTypeString() const {
    switch (_type) {
        case OK: return "OK";
        case MISSED_BALL: return "Missed Ball";
        case INVALID_CUE_PARAMS: return "Invalid Cue Params";
        case POOLFIZ_ERROR: return "Poolfiz Error";
        case UNKNOWN_ERROR: return "Unknown Error";
        default: assert(false);
                 return NULL;
    }
}

ostream& operator<<(ostream &out, const ShotParams &rhs) {
    out << "(a: " << setw(8) << rhs.a;
    out << ", b: " << setw(8) << rhs.b;
    out << ", theta: " << setw(8) << rhs.theta;
    out << ", phi: " << setw(8) << rhs.phi;
    out << ", v: " << setw(8) << rhs.v << ")";
    return out;
}

ostream& operator<<(ostream &out, const Point &rhs) {
    out << "(" << setw(9) << rhs.x << ", " << setw(9) << rhs.y << ")";
    return out;
}

Point Point::rotate(double cos_phi, double sin_phi) const {
    return Point(x*cos_phi - y*sin_phi, 
            x*sin_phi + y*cos_phi);
}

Vector Point::to_v() const {
    return Vector(x, y, 0);
}

Point Point::operator+(const Point& p2) const {
    return Point(x+p2.x, y+p2.y);
}

Point Point::operator-(const Point& p2) const {
    return Point(x-p2.x, y-p2.y);
}

void Point::toStream(ostream &out) const {
    out << setprecision(20) << x << " " << setprecision(20) << y << " ";
}

void Point::fromStream(istream &in) {
    in >> x >> y;;
}

string Point::toString() const {
    ostringstream out;
    toStream(out);
    return out.str();
}

void Point::fromString(const string &s) {
    istringstream in(s);
    fromStream(in);
}

ostream& operator<<(ostream &out, const Vector &rhs) {
    out << "(" << setw(9) << rhs.x << ", " << setw(9) << rhs.y << ", " << setw(1) << rhs.z << ")";
    return out;
}

Vector Vector::operator*(double d) const {
    return Vector(x*d, y*d, z*d);
}


Vector& Pool::Vector::operator *=(double d) 
{
    x*=d; y*=d; z*=d;
    return *this;
}

Vector Vector::operator-(const Vector &fizV) const {
    return Vector(x-fizV.x, y-fizV.y, z-fizV.z);
}

Vector Vector::operator+(const Vector &fizV) const {
    return Vector(x+fizV.x, y+fizV.y, z+fizV.z);
}

Vector Vector::rotateRad(double rad) const {
    return rotate(gsl_sf_cos(rad), gsl_sf_sin(rad));
}

Vector Vector::rotateDeg(double phi) const {
    return rotate(gsl_sf_cos(phi * M_PI/180), gsl_sf_sin(phi * M_PI/180));
}

Vector Vector::rotate(double cos_phi, double sin_phi) const {
    return Vector(x*cos_phi - y*sin_phi, x*sin_phi + y*cos_phi, z);
}

double Vector::dot(const Vector& v2) const {
    return (x*v2.x + y*v2.y + z*v2.z);
}

Vector Vector::cross(const Vector& v2) const {
    return Vector(y*v2.z - z*v2.y, z*v2.x - x*v2.z, x*v2.y - y*v2.x);
}

double Vector::mag() const {
    return sqrt(x*x + y*y + z*z);
}

inline Vector Vector::norm() const {
    //return v_mult(1/v_mag(v), v);
    double m = mag();
#ifdef POOLFIZ_ERRORS    
    if (m == 0) throw "Cannot normalize vector with magnitdue 0";
#endif    
    return Vector(x/m, y/m, z/m);
}

Point Vector::to_p() const {
    return Point(x,y);
}

void Vector::toStream(ostream &out) const {
    out << setprecision(20) << x << " "
        << setprecision(20) << y << " "
        << setprecision(20) << z << " ";
}

void Vector::fromStream(istream &in) {
    in >> x >> y >> z;
}

string Vector::toString() const {
    ostringstream out;
    toStream(out);
    return out.str();
}

void Vector::fromString(const string &s) {
    istringstream in(s);
    fromStream(in);
}

ostream& operator<<(ostream &out, const Ball &rhs) {
    out << "[ID:" << setw(8) << rhs.getIDString() 
        << " " << setw(11) << rhs.getStateString() 
        << " R:" << rhs.getPos() 
        << " V:" << rhs.getVelocity() 
        << " W:" << rhs.getSpin() << "]";
    return out;
}

const Ball& Ball::operator=(const Ball & rhs) {
    if (this != &rhs) {
        radius=rhs.radius;
        r=rhs.r;
        state=rhs.state;
        type=rhs.type;
        v=rhs.v;
        w=rhs.w;
    }
    return *this;
}

string Ball::getIDString() const {
    switch(type) {
        case CUE: return "Cue";
        case ONE: return "One";
        case TWO: return "Two";
        case THREE: return "Three";
        case FOUR: return "Four";
        case FIVE: return "Five";
        case SIX: return "Six";
        case SEVEN: return "Seven";
        case EIGHT: return "Eight";
        case NINE: return "Nine";
        case TEN: return "Ten";
        case ELEVEN: return "Eleven";
        case TWELVE: return "Tweleve";
        case THIRTEEN: return "Thirteen";
        case FOURTEEN: return "Fourteen";
        case FIFTEEN: return "Fifteen";
        case UNKNOWN_ID: return "Unknown ID";
        default: return "defaulted id";
                 assert(false);
                 return NULL;
    }
}

string Ball::getStateString() const {
    switch(state) {
        case NOTINPLAY: return "Not In Play";
        case STATIONARY: return "Stationary";
        case SPINNING: return "Spinning";
        case SLIDING: return "Sliding";
        case ROLLING: return "Rolling";
        case POCKETED_SW: return "Pocketed SW";
        case POCKETED_W: return "Pocketed W";
        case POCKETED_NW: return "Pocketed NW";
        case POCKETED_NE: return "Pocketed NE";
        case POCKETED_E: return "Pocketed E";
        case POCKETED_SE: return "Pocketed SE";
        case SLIDING_SPINNING: return "Sliding-Spinning";
        case ROLLING_SPINNING: return "Rolling-Spinnning";
        case UNKNOWN_STATE: return "Unknown State";
        default: return "defaulted state";
                 assert(false);
                 return NULL;
    }
}

inline bool Ball::overlaps(const Ball &other, double epsilon) const {
    double dx=r.x - other.r.x;
    double dy=r.y - other.r.y;
    //double r=radius+other.radius;
    return /*r*r*/4*BALL_RADIUS*BALL_RADIUS-dx*dx-dy*dy>epsilon;
}

inline double Ball::dist(const Ball &other) const {
    return sqrt( pow((r.x - other.r.x), 2.0) + pow((r.y - other.r.y), 2.0) );
}

void Ball::updateState(bool VERBOSE) {
    //handleStateChange(state.getBall(event->getBall1()));

    if (!isInPlay()) return;
    State newState;
    Vector u = ( v + ((Vector(0, 0, 1).cross(w))*radius) ); //"Relative velocity"
    Vector new_v = v;
    Vector new_w = w;

    if (VERBOSE)
        cerr << "Attempting state change ball " << getIDString() <<"; u is " << u << " (mag = " << u.mag() << "); v is " << v << "; w is " << w << "; old state " << getStateString() <<  endl;

    if(!Utils::vzero(v.mag()))
        //|| !Utils::vzero(u.mag()) )// > Utils::EPSILON*radius)  //The ball is moving
    {
        if(!Utils::vzero(u.mag()) )// > Utils::EPSILON*radius)
            newState = SLIDING;
        else {
            newState = ROLLING;
        }
    }
    else
    {
        if(!Utils::vzero(w.z)) { 
            newState = SPINNING;
            new_w.x = 0;
            new_w.y = 0;
            //new_v = Vector();
        } else {
            newState = STATIONARY;
            new_v = Vector();
            new_w = Vector();
        }
    }

    if (VERBOSE && state != newState) {
        Ball dummy;
        dummy.setState(newState);
        cerr << "Ball " << getIDString() << " transitioned from " << getStateString() << " to " << dummy.getStateString() << endl;
    }
    state = newState;
    v = new_v;
    w = new_w;
}

void Ball::toStream(ostream &out) const {
    out << radius << " " 
        << (int)state << " "
        << (int)type << " ";
    r.toStream(out);
}

void Ball::fromStream(istream &in) {
    in >> radius;
    int s,t;
    in >> s; state = (State)s;
    in >> t; type = (Type)t;
    r.fromStream(in);
    v = Vector();
    w = Vector();
}

string Ball::toString() const {
    ostringstream out;
    toStream(out);
    return out.str();
}

void Ball::fromString(const string &s) {
    istringstream in(s);
    fromStream(in);
}

Table::Table( double length, 
        double width, 
        double cornerPocketWidth, 
        double sidePocketWidth, 
        double muS, 
        double muR, 
        double muSp, 
        double railHeight, 
        double cueLength, 
        double railVelDampingX, 
        double railVelDampingY, 
        double railSpinDamping, 
        double railZSpinAngleAdj, 
        double railVelAngleAdj) :
    _muS(muS), 
    _muR(muR), 
    _muSp(muSp), 
    _railHeight(railHeight), 
    _cueLength(cueLength), 
    _railVelDampingX(railVelDampingX), 
    _railVelDampingY(railVelDampingY), 
    _railSpinDamping(railSpinDamping), 
    _railZSpinAngleAdj(railZSpinAngleAdj), 
    _railVelAngleAdj(railVelAngleAdj) 
{
    calculateTableDimensions( length, width, cornerPocketWidth, sidePocketWidth);
}

Table::Table( const Table &rhs) :
    _muS(rhs._muS),
    _muR(rhs._muR),
    _muSp(rhs._muSp),
    _railHeight(rhs._railHeight),
    _cueLength(rhs._cueLength),
    _railVelDampingX(rhs._railVelDampingX),
    _railVelDampingY(rhs._railVelDampingY),
    _railSpinDamping(rhs._railSpinDamping),
    _railZSpinAngleAdj(rhs._railZSpinAngleAdj),
    _railVelAngleAdj(rhs._railVelAngleAdj),
    _length(rhs._length),
    _width(rhs._width),
    _headString(rhs._headString),
    _footSpot(rhs._footSpot),
    _SWpocketLeft(rhs._SWpocketLeft),
    _SWpocketRight(rhs._SWpocketRight),
    _SWpocketCenter(rhs._SWpocketCenter),
    _WpocketLeft(rhs._WpocketLeft),
    _WpocketRight(rhs._WpocketRight),
    _WpocketCenter(rhs._WpocketCenter),
    _NWpocketLeft(rhs._NWpocketLeft),
    _NWpocketRight(rhs._NWpocketRight),
    _NWpocketCenter(rhs._NWpocketCenter),
    _NEpocketLeft(rhs._NEpocketLeft),
    _NEpocketRight(rhs._NEpocketRight),
    _NEpocketCenter(rhs._NEpocketCenter),
    _EpocketLeft(rhs._EpocketLeft),
    _EpocketRight(rhs._EpocketRight),
    _EpocketCenter(rhs._EpocketCenter),
    _SEpocketLeft(rhs._SEpocketLeft),
    _SEpocketRight(rhs._SEpocketRight),
    _SEpocketCenter(rhs._SEpocketCenter) {}

    void Table::calculateTableDimensions( double length, double width, double cornerPocketWidth, double sidePocketWidth ) {
        _length = length;
        _width = width;
        double diagonalWidth = cornerPocketWidth / sqrt(2.0);
        double straightWidth = sidePocketWidth / 2.0;
        _SWpocketLeft = Point( diagonalWidth, 0.0 );
        _SWpocketRight = Point( 0.0, diagonalWidth );
        _SWpocketCenter = Point( diagonalWidth/2.0, diagonalWidth/2.0 );
        _WpocketLeft = Point( 0.0, length/2.0 - straightWidth );
        _WpocketRight = Point( 0.0, length/2.0 + straightWidth );
        _WpocketCenter = Point( 0.0, length/2.0 );
        _NWpocketLeft = Point( 0.0, length - diagonalWidth );
        _NWpocketRight = Point( diagonalWidth, length );
        _NWpocketCenter = Point( diagonalWidth/2.0, length - diagonalWidth/2.0 );
        _NEpocketLeft = Point( width - diagonalWidth, length );
        _NEpocketRight = Point( width, length - diagonalWidth );
        _NEpocketCenter = Point( width - diagonalWidth/2.0, length - diagonalWidth/2.0 );
        _EpocketLeft = Point( width, length/2.0 + straightWidth );
        _EpocketRight = Point( width, length/2.0 - straightWidth );
        _EpocketCenter = Point( width, length/2.0 );
        _SEpocketLeft = Point( width, diagonalWidth );
        _SEpocketRight = Point( width - diagonalWidth, 0.0 );
        _SEpocketCenter = Point( width - diagonalWidth/2.0, diagonalWidth/2.0 );		

        _headString = length * 3/4;
        _footSpot = Point(width/2, length/4);
    }

const Table & Table::operator=(const Table & rhs) {
    if (this != &rhs) {
        _muS = rhs._muS;
        _muR = rhs._muR;
        _muSp = rhs._muSp;
        _railHeight = rhs._railHeight;
        _cueLength = rhs._cueLength;
        _railVelDampingX = rhs._railVelDampingX;
        _railVelDampingY = rhs._railVelDampingY;
        _railSpinDamping = rhs._railSpinDamping;
        _railZSpinAngleAdj = rhs._railZSpinAngleAdj;
        _railVelAngleAdj = rhs._railVelAngleAdj;
        _length = rhs._length;
        _width = rhs._width;
        _headString = rhs._headString;
        _footSpot = rhs._footSpot;
        _SWpocketLeft = rhs._SWpocketLeft;
        _SWpocketRight = rhs._SWpocketRight;
        _SWpocketCenter = rhs._SWpocketCenter;
        _WpocketLeft = rhs._WpocketLeft;
        _WpocketRight = rhs._WpocketRight;
        _WpocketCenter = rhs._WpocketCenter;
        _NWpocketLeft = rhs._NWpocketLeft;
        _NWpocketRight = rhs._NWpocketRight;
        _NWpocketCenter = rhs._NWpocketCenter;
        _NEpocketLeft = rhs._NEpocketLeft;
        _NEpocketRight = rhs._NEpocketRight;
        _NEpocketCenter = rhs._NEpocketCenter;
        _EpocketLeft = rhs._EpocketLeft;
        _EpocketRight = rhs._EpocketRight;
        _EpocketCenter = rhs._EpocketCenter;
        _SEpocketLeft = rhs._SEpocketLeft;
        _SEpocketRight = rhs._SEpocketRight;
        _SEpocketCenter = rhs._SEpocketCenter;
    }
    return *this;
}

const Point & Table::getPocketCenter(Pocket pocket) const {
    switch(pocket) {
        case SW: return _SWpocketCenter;
        case W: return _WpocketCenter;
        case NW: return _NWpocketCenter;
        case NE: return _NEpocketCenter;
        case E: return _EpocketCenter;
        case SE: return _SEpocketCenter;
        default: assert(false);
                 Point* dummy = new Point();
                 return *dummy;
    }
}

const Point & Table::getPocketRight(Pocket pocket) const {
    switch(pocket) {
        case SW: return _SWpocketRight;
        case W: return _WpocketRight;
        case NW: return _NWpocketRight;
        case NE: return _NEpocketRight;
        case E: return _EpocketRight;
        case SE: return _SEpocketRight;
        default: assert(false);
                 Point* dummy = new Point();
                 return *dummy; 
    }
}

const Point & Table::getPocketLeft(Pocket pocket) const {
    Point left = Point();
    switch(pocket) {
        case SW: return _SWpocketLeft;
        case W: return _WpocketLeft;
        case NW: return _NWpocketLeft;
        case NE: return _NEpocketLeft;
        case E: return _EpocketLeft;
        case SE: return _SEpocketLeft;
        default: assert(false);
                 Point* dummy = new Point();
                 return *dummy;
    }
}

Ball::State Table::stateFromPocket(Pocket pocket) {
    switch(pocket)
    {
        case SW: return Ball::POCKETED_SW;
        case W: return Ball::POCKETED_W;
        case NW: return Ball::POCKETED_NW;
        case NE: return Ball::POCKETED_NE;
        case E: return Ball::POCKETED_E;
        case SE: return Ball::POCKETED_SE;
        default:
                 return Ball::UNKNOWN_STATE;
    }
}

Table::Pocket Table::pocketFromBndId(BoundaryId bnd) {
    switch(bnd)
    {
        case SW_POCKET: return SW;
        case W_POCKET: return W;
        case NW_POCKET: return NW;
        case NE_POCKET: return NE;
        case E_POCKET: return E;
        case SE_POCKET: return SE;
        default:
                        return UNKNOWN_POCKET;
    }
}

Table::BoundaryId Table::bndIdFromPocket(Pocket pocket) {
    switch(pocket)
    {
        case SW: return SW_POCKET;
        case W: return W_POCKET;
        case NW: return NW_POCKET;
        case NE: return NE_POCKET;
        case E: return E_POCKET;
        case SE: return SE_POCKET;
        default:
                 return UNKNOWN_BOUNDARY;
    }
}

string Table::boundaryName(BoundaryId boundary) {
    switch(boundary)
    {
        case SW_POCKET:return "SW Pocket";
        case SW_RAIL:return "SW Rail";
        case W_POCKET:return "W Pocket";
        case NW_RAIL:return "NW Rail";
        case NW_POCKET:return "NW Pocket";
        case N_RAIL:return "N Rail";
        case NE_POCKET:return "NE Pocket";
        case NE_RAIL:return "NE Rail";
        case E_POCKET:return "E Pocket";
        case SE_RAIL:return "SE Rail";
        case SE_POCKET:return "SE Pocket";
        case S_RAIL: return "S Rail";
        case UNKNOWN_BOUNDARY:return "Unknown Boundary";
        default: assert(false);
                 return NULL;
    }
}

string Table::pocketName(Pocket pocket) {
    switch(pocket)
    {
        case SW:return "SW Pocket";
        case W:return "W Pocket";
        case NW:return "NW Pocket";
        case NE:return "NE Pocket";
        case E:return "E Pocket";
        case SE:return "SE Pocket";
        case UNKNOWN_POCKET:return "Unknown Pocket";
        default: assert(false);
                 return NULL;
    }
}

bool Event::operator<(const Event &other) const {
    return _time < other._time;
}

bool Event::eventCmp(const Event* event1, const Event* event2) {
    return *event1 < *event2;
}

bool Event::relatedTo(const Event &other) const {
    return other.involvesBall(_ball1);
}

bool Event::involvesBall(Ball::Type b) const {
    return (_ball1 == b);
}

bool BallCollisionEvent::relatedTo(const Event& other) const {
    return other.involvesBall(_ball1) || other.involvesBall(_ball2);
}

bool BallCollisionEvent::involvesBall(Ball::Type b) const {
    return ((_ball1 == b) || (_ball2 == b));
}

void StateChangeEvent::doHandle(TableState &ts, bool VERBOSE) const {
    ts.getBall(_ball1).updateState(VERBOSE);
}

void BallCollisionEvent::doHandle(TableState &ts, bool VERBOSE) const {
    Ball &ball1 = ts.getBall(_ball1);
    Ball &ball2 = ts.getBall(_ball2);

    //TODO: why is this an error?
    //if(ball2.getID() == Ball::CUE)
    //cerr << "Error in Shot::handleEvent! Cue ball chosen instead of other ball" << endl;

    //handleBallCollision(state.getBall(ballEvent->getBall1()), state.getBall(ballEvent->getBall2()));

    if(VERBOSE)
    {
        cerr << "Handling ball/ball collision" << endl;

        cerr << "Ball 1: " << ball1.getIDString() << ", " << ball1.getVelocity() << ", "
            << ball1.getSpin() << endl;
        cerr << "Ball 2: " << ball2.getIDString() << ", " << ball2.getVelocity() << ", "
            << ball2.getSpin() << endl;
    }

    if (ball1.getState() == Ball::STATIONARY) {
        if (VERBOSE) {
            cerr << "Collision by a stationary ball?" << endl;
            return;
        }
    }

    //First, align the two balls to be on the x-axis.
    double phi = Utils::angle(ball1.getPos(), ball2.getPos(), Utils::START_ZERO, Utils::DEGREES);
    ball1.setVelocity(ball1.getVelocity().rotateDeg(-phi));
    ball2.setVelocity(ball2.getVelocity().rotateDeg(-phi));
    ball1.setSpin(ball1.getSpin().rotateDeg(-phi));
    ball2.setSpin(ball2.getSpin().rotateDeg(-phi));


    double mu_B = 0.01;
    double v0x = ball1.getVelocity().x;// - ball2.getVelocity().x;

    Vector vBB; //relative velocity of balls at contact point
    vBB = (ball1.getVelocity() + (ball2.getVelocity()*(-1)) + (Vector(1, 0, 0).cross(ball1.getSpin() + ball2.getSpin())*-0.028575) );
    double vBB_mag = vBB.mag();

    Vector new_obj_v; //referring to ball2
    new_obj_v.y = ball2.getVelocity().y + mu_B*v0x*vBB.y/vBB_mag;
    new_obj_v.x = ball1.getVelocity().x;
    new_obj_v.z = 0;

    Vector new_cue_v; //ball1
    new_cue_v.y = ball1.getVelocity().y - mu_B*v0x*vBB.y/vBB_mag;
    new_cue_v.x = ball2.getVelocity().x;
    new_cue_v.z = 0;

    Vector new_obj_w, new_cue_w;
    //TODO: Test behavior of spinning, rolling, sliding/rolling cases
    if(ball2.getState() == Ball::STATIONARY || ball2.getState() == Ball::SPINNING)
    {
        //Follow Marlow...
        Vector w_to_add = ((Vector(1, 0, 0).cross(vBB))*(-2.5*mu_B*v0x/vBB_mag/ball1.getRadius()/1.000575));

        new_obj_w = ( ball2.getSpin() + w_to_add );
        new_cue_w = ( ball1.getSpin() + (w_to_add*(-1)) );
    }
    else
    {
        //Apparently Poolfiz just makes the balls roll...
        new_obj_w.x = -new_obj_v.y/ball1.getRadius();
        new_obj_w.y = new_obj_v.x/ball1.getRadius();
        new_obj_w.z = ball2.getSpin().z;

        new_cue_w.x = -new_cue_v.y/ball1.getRadius();
        new_cue_w.y = new_cue_v.x/ball1.getRadius();
        new_cue_w.z = ball1.getSpin().z;
    }

    ball1.setVelocity(new_cue_v);
    ball2.setVelocity(new_obj_v);
    ball1.setSpin(new_cue_w);
    ball2.setSpin(new_obj_w);

    //Finally, unrotate the balls to their original collision angle.
    ball1.setVelocity(ball1.getVelocity().rotateDeg(phi));
    ball2.setVelocity(ball2.getVelocity().rotateDeg(phi));
    ball1.setSpin(ball1.getSpin().rotateDeg(phi));
    ball2.setSpin(ball2.getSpin().rotateDeg(phi));

    //ball1.checkEpsilons();
    //ball2.checkEpsilons();
    if(Utils::vzero(ball1.getVelocity().mag()) && Utils::vzero(ball1.getSpin().mag()))
    {
        ball1.setVelocity(Vector());
        ball1.setSpin(Vector(0, 0, ball1.getSpin().z));
    }
    if(Utils::vzero(ball2.getVelocity().mag()) && Utils::vzero(ball2.getSpin().mag()))
    {
        ball2.setVelocity(Vector());
        ball2.setSpin(Vector(0, 0, ball2.getSpin().z));
    }

    if(VERBOSE)
    {
        cerr << "Ball 1: " << ball1.getVelocity() << ", "
            << ball1.getSpin() << endl;
        cerr << "Ball 2: " << ball2.getVelocity() << ", "
            << ball2.getSpin() << endl;
    }

    //handleStateChange(ball1);
    //handleStateChange(ball2);
    ts.getBall(_ball1).updateState(VERBOSE);
    ts.getBall(_ball2).updateState(VERBOSE);
}

void RailCollisionEvent::doHandle(TableState &ts, bool VERBOSE) const {
    Ball &ball = ts.getBall(_ball1);

    if(VERBOSE)
        cerr << "Handling rail collision of " << ball.getIDString() 
            << " with " << Table::boundaryName(_rail) << endl;
    //handleRailCollision(state.getBall(railEvent->getBall1()), railEvent->getRail());

    Vector v = ball.getVelocity();
    Vector newVelocity;

    if(_rail == Table::N_RAIL || _rail == Table::S_RAIL) //Horizontal rails
        newVelocity = Vector(0.9*v.x, -0.9*v.y, 0);
    else
        newVelocity = Vector(-0.9*v.x, 0.9*v.y, 0);
    ball.setVelocity(newVelocity);

    //Spin is just divided by 10
    ball.setSpin(ball.getSpin()*0.1);

    //The ball will usually be sliding at this point, whatever it was doing before
    //handleStateChange(ball);
    ts.getBall(_ball1).updateState(VERBOSE);

    if(VERBOSE)
    {
        cerr << ball.getVelocity() << ", "
            << ball.getSpin() << endl;
    }
}

void PocketedEvent::doHandle(TableState &ts, bool VERBOSE) const {
    Ball &ball = ts.getBall(_ball1);

    ball.setState(Table::stateFromPocket(_pocket));
    if(VERBOSE)
        cerr << "Ball " << ball.getIDString() << " pocketed in "
            << Table::pocketName(_pocket) << endl;

    ball.setVelocity(Vector());
    ball.setSpin(Vector());
}

void CueStrikeEvent::doHandle(TableState &ts, bool VERBOSE) const {
    Ball &ball = ts.getBall(_ball1);

    double m = Ball::BALL_MASS;
    double M = Table::CUE_MASS;
    double cos = gsl_sf_cos(_params.theta * M_PI/180);
    double sin = gsl_sf_sin(_params.theta * M_PI/180);
    double c = fabs(sqrt(ball.getRadius()*ball.getRadius() - _params.a*_params.a - _params.b*_params.b));
    //TODO: Check that c is sane

    //correctional constant to match poolfiz
    double c1 = 0.999425;

    //F is the force of the cue stick on the ball
    //Equation (1) in Leckie and Greenspan
    double F = (2*m*_params.v) / (1 + (m/M) + c1*(5.0/(2.0*ball.getRadius()*ball.getRadius()))*(_params.a*_params.a + _params.b*_params.b*cos*cos + c*c*sin*sin - 2*_params.b*c*cos*sin));
    if(VERBOSE)
        printf("F = %g\n", F);
    //Now, set the velocity of the cue ball
    //In i,j,k coordinate system:
    //velocity is (0, -F*cos/m, 0), making no-jumps assumption (2)
    //angular velocity is (1/I)*(-c*F*sin + b*F*cos, a*F*sin, -a*F*cos) (3)

    //The following equations ((2) and (3) from Leckie/Greenspan) are in a
    //ball-centric coordinate system (ijk). j is parallel to the projection of
    //the cue onto the table; k is the same as z; and the ball's center is at
    //(0, 0, 0). v_rotate returns this orientation to that of the table.
    ball.setVelocity(Vector(0, -F*cos/m, 0).rotateDeg(_params.phi + 90.0));
    ball.setSpin((Vector(-c*F*sin + _params.b*F*cos, _params.a*F*sin, -_params.a*F*cos)*(1/Table::I)).rotateDeg(_params.phi + 90.0));

    //Ball is most likely sliding now
    //handleStateChange(ball);
    ts.getBall(_ball1).updateState(VERBOSE);

    if(VERBOSE)
    {
        printf("Just struck the cue ball; cue ball has state:\n");
        cerr << ball.getVelocity() << ", "
            << ball.getSpin() << endl;

        //cerr << state.getBall(CUE) << endl;
    }
}

void Event::copyBalls(TableState &ts) {
    if (_ball1Data != NULL) delete _ball1Data;
    _ball1Data = new Ball(ts.getBall(_ball1));
}

void BallCollisionEvent::copyBalls(TableState &ts) {
    Event::copyBalls(ts);
    if (_ball2Data != NULL) delete _ball2Data;
    _ball2Data = new Ball(ts.getBall(_ball2));
}

string Event::toString() const {
    ostringstream out;
    dump(out);
    return out.str();
}

ostream& operator<<(ostream &out, const Event &rhs) {
    out << "<";
    rhs.dump(out);
    out << ">";
    return out;
}

ostream& Event::dump(ostream &out) const {
    out << "T:" << setw(8) << getTime() << " Type:" << setw(14) << getTypeString() << " ";
    if (_ball1Data != NULL) 
        out << *_ball1Data;
    else {
        Ball dummy(_ball1);
        out << dummy.getIDString();
    }
    return out;
}

ostream& BallCollisionEvent::dump(ostream &out) const {
    Event::dump(out);
    out << " ";
    if (_ball2Data != NULL) 
        out << *_ball2Data;
    else {
        Ball dummy(_ball2);
        out << dummy.getIDString();
    }
    return out;
}

ostream& RailCollisionEvent::dump(ostream &out) const {
    Event::dump(out);
    out << " Rail:" << setw(7) << Table::boundaryName(getRail());
    return out;
}

ostream& PocketedEvent::dump(ostream &out) const {
    Event::dump(out);
    out << " Pocket:" << setw(7) << Table::pocketName(getPocket());
    return out;
}

ostream& CueStrikeEvent::dump(ostream &out) const {
    Event::dump(out);
    out << " Params:" << getParams();
    return out;
}

ostream& operator<<(ostream &out, Shot &rhs) {
    out << "Events {" << endl;
    for (unsigned int i=0; i<rhs.getEventList().size(); i++) {
        if (i!=0) {
            out << " dx: " << (rhs.getEventList()[i]->getTime()-rhs.getEventList()[i-1]->getTime()) << endl;
        }
        out << " " << *(rhs.getEventList()[i]) << endl;
    }
    out <<"}";
    return out;
}

Shot::~Shot() {
    for (vector<Event*>::iterator iter = shotEvents.begin(); iter != shotEvents.end(); ++iter)
        delete *iter;
}

double Shot::getDuration() const {
    return shotEvents.back()->getTime();
}

void Shot::simulateShot(TableState &state, const ShotParams &sp, bool fullSim, bool verbose, bool errors)
{
    /*
#ifndef NDEBUG
dbgState=state;
dbgParams.a=a;
dbgParams.b=b;
dbgParams.theta=theta;
dbgParams.phi=phi;
dbgParams.V=V;
signal(SIGABRT,&error);
#endif
*/
    VirtualStopwatch::step();
    //int ballevents[50]={0};
    long eventcounter=5000;
    double a = sp.a;
    double b = sp.b;
    if (verbose) {
        cerr<< "Executing shot " << sp << " on table " << state<< endl;
    }
    gsl_set_error_handler_off();
    a /= -1000.0;
    b /= 1000.0;
    VERBOSE = verbose;
    CHECK_ERRORS = errors;

    //Put the table state into a form we can use
    //TableState *newState = &state;

    double curTime = 0;
    list<Event*> futureEvents; //Events that may or may not happen after curTime

    //If the shot isn't physically possible, stop.
    if(TableState::OK_PRECONDITION != state.isPhysicallyPossible(sp,VERBOSE))
    {
        BadShotException e(BadShotException::INVALID_CUE_PARAMS);
        if (VERBOSE) {
            cerr << "Fastfiz exited early - invalid cue params " << state.isPhysicallyPossible(sp) << endl;
        }
        throw e;
    }
    //TODO: Signal a miscue when appropriate

    string s = state.toString();

    //Add the cue strike to the list of future events
    futureEvents.push_back(new CueStrikeEvent(ShotParams(a, b, sp.theta, sp.phi, sp.v)));

    while(!futureEvents.empty() && (--eventcounter))
    { 
        VirtualStopwatch2::step();
        //Figure out which event comes first
        list<Event*>::iterator currentEventItr = min_element(futureEvents.begin(), futureEvents.end(), Event::eventCmp);
        Event *currentEvent = *currentEventItr;
        futureEvents.erase(currentEventItr);

        if (VERBOSE) {
            list<Event*>::iterator nextEventItr = min_element(futureEvents.begin(), futureEvents.end(), Event::eventCmp);
            cerr << "Counter: " << eventcounter << endl;
            cerr << "THIS EVENT " << *currentEvent << endl;

            if (nextEventItr != futureEvents.end())
            {
                Event *nextEvent = *nextEventItr;
                if (nextEvent) {
                    cerr << "NEXT EVENT " << *nextEvent <<endl;
                    double delta=(nextEvent->getTime() - currentEvent->getTime());
                    cerr << "Time delta: " << delta << endl;
                    if (Utils::fequal(delta,0)) {
                        cerr << "VERY LOW DELTA!" << endl;
                    }
                }
            }
        }

        /*
           if (eventcounter<3000) {
           state.addNoise(1e-7);
           if (VERBOSE) {
           cerr << "INFINITE LOOP AVOIDANCE -- adding random noise!" << endl;
           }
           }*/

        //cerr << "future events:" << futureEvents.size() << endl;
        //cerr << *currentEvent << endl;
        //cerr << "FUTURE LOOP\n";
        //currentEvent->dump();

        //Sanity check
        while(!(currentEvent->getTime() >= curTime))
        {
            if (VERBOSE) {
                cerr << "sanity check triggered" << endl;
                cerr << "cur time: " << setprecision(30) << curTime << endl;
                cerr << "evt time: " << setprecision(30) << currentEvent->getTime() << endl;
                cerr << *currentEvent << endl;
            }
            POOLFIZ_ERROR();
            delete currentEvent; currentEvent=NULL;
            currentEventItr = min_element(futureEvents.begin(), futureEvents.end(), Event::eventCmp);
            if(currentEventItr == futureEvents.end())
                return;
            currentEvent = *currentEventItr;
            futureEvents.erase(currentEventItr);
        }

        //cerr << "FUTURE LOOP after selection\n";
        //currentEvent->dump();

        //Update each ball's position to the next event's time
        if (VERBOSE) {
            cerr << "before update" << state << endl; 
        }


        updateTime(state, curTime, currentEvent->getTime());

        /*
           if (VERBOSE) {
           cerr << "is valid: " << state.isValidBallPlacement(VERBOSE) << endl;
           cerr << s << endl;
           }
           */

        if (VERBOSE) {
            //cerr << "be [" << *currentEvent << "] = " << ballevents[currentEvent->getBall1()]++ << endl;
            //for (int bt=0; bt<16; bt++) {
            //  cerr << ballevents[bt] << ";";
            //}
            cerr << endl;
            cerr << "after update and before handle" << state << endl;
            cout <<endl;
            state.toStream(cout);
            cout <<endl;
            printOverlap(state);
            if (state.isValidBallPlacement(true) != TableState::OK_PRECONDITION) {
                cerr << "is valid: " << state.isValidBallPlacement(true) << endl;
                cerr << *currentEvent<<endl;
                for (list<Event*>::iterator itr = futureEvents.begin(); itr!=futureEvents.end(); ++itr) {
                    if ((*itr)->getType() == Event::BALL_COLLISION)
                        cerr << *(*itr) << endl;
                }
            } 
        }

#ifdef POOLFIZ_ERRORS
        if (CHECK_ERRORS && state.isValidBallPlacement() != TableState::OK_PRECONDITION) {
            POOLFIZ_ERROR();
        }
#endif

        if (VERBOSE) {
            cerr << "difference in time" << (currentEvent->getTime() - curTime) << endl;
        }
        curTime = currentEvent->getTime();


        //If we're just looking for the first ball hit, check for that...
        if(!fullSim && currentEvent->getType() == Event::BALL_COLLISION)
        {
            shotEvents.push_back(currentEvent);
            //delete currentEvent; 
            currentEvent=NULL;
            //Free memory before returning
            for(list<Event*>::iterator itr = futureEvents.begin(); itr != futureEvents.end(); ++itr)
            {
                delete *itr;
            }
            return;
        }

        //Now we change the table state according to the event itself
        currentEvent->handle(state, VERBOSE);
        /*
           if (state.isValidBallPlacement() != TableState::OK_PRECONDITION) {
           if (VERBOSE) {
           cerr << "is valid: " << state.isValidBallPlacement(VERBOSE) << endl;
           cerr << s << endl;
           }
           }
           */
#ifdef POOLFIZ_ERRORS
        if (CHECK_ERRORS && state.isValidBallPlacement() != TableState::OK_PRECONDITION) {
            POOLFIZ_ERROR();
        }
#endif
        if (VERBOSE) {
            cerr << "after handle" << state << endl;
            if (state.isValidBallPlacement(true) != TableState::OK_PRECONDITION) {
                cerr << "is valid: " << state.isValidBallPlacement() << endl;
                printOverlap(state);
            }
        }
        /*if (ballevents[currentEvent->getBall1()]++>200) {
        // STOP BALL
        Ball &b=state.getBall(currentEvent->getBall1());
        if (VERBOSE) {
        cerr << "Timeout failsaife triggered for ball " << b << ". Stopping ball." << endl;
        }
        b.setVelocity(Vector(0,0,0));
        b.setSpin(Vector(0,0,0));
        b.setState(Ball::STATIONARY);
        }*/
        //handleEvent(state, currentEvent);

        /*
           if (LOGFILE) {
           char comment[100];
           sprintf(comment,"Curr time = %lf",curTime);
           TableUtils::writeLogComment(LOGFILE,comment);
           TableUtils::writeLog(LOGFILE,state);
           }
           */

        //Check for new events ("simple" method: wipe futureEvents and check everything)
        /*for(list<Event*>::iterator itr = futureEvents.begin(); itr != futureEvents.end(); itr++)
          delete *itr;
          futureEvents.clear();
          addFutureEvents(table, newState, curTime, futureEvents);*/

        //cerr << "FUTURE LOOP after handle\n";
        //currentEvent->dump();

        //Check for new events (smarter method: update futureEvents only for those
        //balls involved in the last event)
        if (VERBOSE) {
            //cerr << "before remove" << endl;
            //printEvents(futureEvents);
        }
#ifdef ALL_EVENTS
        removeRelatedEvents(futureEvents, NULL); //REMOVE ALL EVENTS
        addRelatedFutureEvents(state,curTime,futureEvents,Ball::UNKNOWN_ID,Ball::UNKNOWN_ID); // ADD ALL EVENTS
#else        
        removeRelatedEvents(futureEvents, currentEvent);
        if (VERBOSE) {
            //cerr << "after remove and before addRelatedFuture" << endl;
            //printEvents(futureEvents);
        }
        //cerr << "FUTURE LOOP after remove\n";
        //currentEvent->dump();

        addRelatedFutureEvents(state, curTime, futureEvents, currentEvent->getBall1(), currentEvent->getBall2());
#endif        
        if (VERBOSE) {
            //cerr << "after addRelatedFuture" << endl;
            //printEvents(futureEvents);
        }

        shotEvents.push_back(currentEvent);
        //delete currentEvent; 
        currentEvent=NULL;
    }
    if (!eventcounter) {
        BadShotException e(BadShotException::POOLFIZ_ERROR);
        //if (VERBOSE) {
        cerr << "Fastfiz exited early - infinite loop. " << endl;
        cerr << "sp= " << sp << endl;
        //}
        throw e;
    }
}

void Shot::updateTime(TableState &state, double oldTime, double newTime)
{
    if(VERBOSE)
        cerr << "Time " << newTime << ": ";
    if(oldTime > newTime)
        throw "Error! Trying to step to a previous time"; 
    //printf("Error! Trying to step from time %g to time %g\n", oldTime, newTime);
    for(vector<Ball>::iterator i = state.getBegin(); i != state.getEnd(); ++i)
    {
        updateBall(state.getTable(), *i, oldTime, newTime);
    }
    state.fixOverlap(VERBOSE);
}

void Shot::updateBall(const Table &table, Ball &ball, double oldTime, double newTime)
{
    if(!(ball.getState() == Ball::SPINNING || ball.getState() == Ball::ROLLING || ball.getState() == Ball::SLIDING))
        return;

    if((ball.getID() == Ball::CUE || ball.getID() == Ball::ONE) && VERBOSE)
    {
        cerr << "Ball " << ball.getID() << " is stepping forward from time " << oldTime << " to time " << newTime << "." << endl;
        cerr << "At " << oldTime << ": " << ball << endl;
    }

    Point r = ball.getPos();
    Vector v = ball.getVelocity();
    Vector w = ball.getSpin();
    double t = newTime - oldTime; // LOSS OF PRECISION!
    double mu_sp = table.getMuSpinning();
    double mu_r = table.getMuRolling();
    double mu_s = table.getMuSliding();
    Vector u_o;
    Vector vNorm;
    Point new_r = r;
    Vector new_v = v;
    Vector new_w = w;

    double cos_phi = 1.0, sin_phi = 0.0;
    if(v.mag() > 0)
    {
        vNorm = v.norm();
        //Used for shifting in and out of a ball-centric reference frame
        cos_phi = vNorm.x;
        sin_phi = vNorm.y;

        if((ball.getID() == Ball::CUE || ball.getID() == Ball::ONE) && VERBOSE) {
            cerr << " before rotate v" << v << " w" << w << endl;
        }
        //Rotate by -phi; v should end up along the x-axis, so the equations will work
        v = v.rotate(cos_phi, -sin_phi);
        w = w.rotate(cos_phi, -sin_phi);
        if((ball.getID() == Ball::CUE || ball.getID() == Ball::ONE) && VERBOSE) {
            cerr << " after rotate v" << v << " w" << w << endl;
        }
    }

    //New velocity, spin, etc. depends on whether it's sliding, rolling, or spinning... or still
    switch(ball.getState())
    {
        case Ball::SPINNING:
            new_w.x = 0;
            new_w.y = 0;
            new_w.z = updateSpinning(w.z, t, mu_sp, ball.getRadius(), false);
            break;

        case Ball::SLIDING:
            //Equations (4)-(8) in Leckie/Greenspan 
            if((ball.getID() == Ball::CUE || ball.getID() == Ball::ONE) && VERBOSE) {
                cerr << "current ball v" << v << endl;
                cerr << "current ball w" << w << endl;
            }

            mu_s = table.getMuSliding(); //should be 0.2
            u_o = (v + (Vector(0,0,1).cross(w)*ball.getRadius())).norm();
            new_r.x = v.mag()*t - (0.5)*mu_s*Table::g*t*t*u_o.x;
            new_r.y = -(0.5)*mu_s*Table::g*t*t*u_o.y;
            new_v = ( v + (u_o*(-mu_s*Table::g*t)) );
            new_w = ( w + (u_o.cross(Vector(0, 0, 1))*(-(5.0*mu_s*Table::g)/(2.0*ball.getRadius())*t)) );
            new_w.z = updateSpinning(w.z, t, mu_sp, ball.getRadius(), true);

            if((ball.getID() == Ball::CUE || ball.getID() == Ball::ONE) && VERBOSE) {
                cerr << "new_v" << new_v << endl;          
                cerr << "new_w" << new_w << endl;
                cerr << "u_o" << u_o << endl;
            }

            new_r = new_r.rotate(cos_phi, sin_phi);
            new_r = r + new_r;
            new_v = new_v.rotate(cos_phi, sin_phi);
            new_w = new_w.rotate(cos_phi, sin_phi);

            break;

        case Ball::ROLLING:
            //Equations (12)-(14) and (8) in Leckie/Greenspan
            vNorm = v.norm();

            mu_r = table.getMuRolling(); //should be 0.015
            new_r = ( (v*t) + (vNorm*(-(0.5)*mu_r*Table::g*t*t)) ).to_p();
            new_v = ( v + (vNorm*(-mu_r*Table::g*t)) );
            new_w = ( w * (new_v.mag()/v.mag()) );
            new_w.z = updateSpinning(w.z, t, mu_sp, ball.getRadius(), false);

            new_r = new_r.rotate(cos_phi, sin_phi);
            new_r = (r + new_r);
            new_v = new_v.rotate(cos_phi, sin_phi);
            new_w = new_w.rotate(cos_phi, sin_phi);
            break;

        default:
            //Not moving or not in play; do nothing
            return;
    }
    if(Utils::fequal(new_v.x, 0))
        new_v.x = 0;
    if(Utils::fequal(new_v.y, 0))
        new_v.y = 0;
    if(Utils::fequal(new_v.z, 0))
        new_v.z = 0;
    if(Utils::fequal(new_w.x, 0))
        new_w.x = 0;
    if(Utils::fequal(new_w.y, 0))
        new_w.y = 0;
    if(Utils::fequal(new_w.z, 0))
        new_w.z = 0;

    ball.setPos(new_r);
    ball.setVelocity(new_v);
    ball.setSpin(new_w);
    ball.updateState(VERBOSE);
    //if(VERBOSE)
    //cerr << "At " << newTime << ": " << ball << endl;
}

double Shot::updateSpinning(double w_z, double t, double mu_sp, double R, bool isSliding)
{
    if(Utils::fequal(w_z, 0))
        return 0;
    //The sliding factor is to match Poolfiz
    double new_w_z = w_z - 5*mu_sp*Table::g*t/(2*Ball::BALL_RADIUS) * (w_z > 0 ? 1 : -1) * (isSliding ? 0.25 : 1);
    if ((new_w_z * w_z) <= 0.0) //Stop it at zero
        new_w_z = 0;
    return new_w_z;
}

/*
   void Shot::handleEvent(TableState &state, Event *event)
   {
   BallCollisionEvent* ballEvent;
   RailCollisionEvent* railEvent;
   PocketedEvent* pocketEvent;
   CueStrikeEvent* cueEvent;
   switch(event->getType())
   {
   case Event::STATE_CHANGE:
   if(VERBOSE)
   cerr << "Handling state change of " << event->getBall1() << endl;
   handleStateChange(state.getBall(event->getBall1()));
   break;

   case Event::BALL_COLLISION:
   ballEvent = dynamic_cast<BallCollisionEvent*>(event);
   if(ballEvent->getBall2() == Ball::CUE)
   cerr << "Error in Shot::handleEvent! Cue ball chosen instead of other ball" << endl;
   if(VERBOSE)
   cerr << "Collision of " << ballEvent->getBall1() << " and " << ballEvent->getBall2() << ": " << endl;
   handleBallCollision(state.getBall(ballEvent->getBall1()), state.getBall(ballEvent->getBall2()));
   break;

   case Event::RAIL_COLLISION:
   railEvent = dynamic_cast<RailCollisionEvent*>(event);
   if(VERBOSE)
   cerr << "Handling rail collision of " << railEvent->getBall1() << endl;
   handleRailCollision(state.getBall(railEvent->getBall1()), railEvent->getRail());
   if(VERBOSE)
   {
   Ball ball1 = state.getBall(railEvent->getBall1());
   cerr << "(" << ball1.getVelocity().x << ", "
   << ball1.getVelocity().y << ", "
   << ball1.getVelocity().z << "), ("
   << ball1.getSpin().x << ", "
   << ball1.getSpin().y << ", "
   << ball1.getSpin().z << ")" << endl;
   }
   break;

   case Event::POCKETED:
   pocketEvent = dynamic_cast<PocketedEvent*>(event);
   handlePocketed(state.getBall(pocketEvent->getBall1()), pocketEvent->getPocket());
   break;

   case Event::CUE_STRIKE:
   cueEvent = dynamic_cast<CueStrikeEvent*>(event);
   handleCueStrike(state.getBall(Ball::CUE), cueEvent->getParams());
   if(VERBOSE)
   {
   printf("Just struck the cue ball; cue ball has state:\n");
   Ball ball1 = state.getBall(Ball::CUE);
   cerr << "(" << ball1.getVelocity().x << ", "
   << ball1.getVelocity().y << ", "
   << ball1.getVelocity().z << "), ("
   << ball1.getSpin().x << ", "
   << ball1.getSpin().y << ", "
   << ball1.getSpin().z << ")" << endl;

//cerr << state.getBall(CUE) << endl;
}
break;

case Event::UNKNOWN_EVENT:
default:
printf("Non-event \"%d\" was received by Shot::handleEvent().\n", event->getType());
break;
}

}

//Given a ball, updates the ball's state to match its velocity and spin.
void Shot::handleStateChange(Ball &ball)
{
    //Figure out the new state
    if (!ball.isInPlay()) return;
    Ball::State newState;
    Vector u = ( ball.getVelocity() + ((Vector(0, 0, 1).cross(ball.getSpin()))*Ball::BALL_RADIUS) ); //"Relative velocity"

    if (VERBOSE)
        cerr << "Attempting state change ball " << ball.getID() <<"; u is " << u << " (mag = " << u.mag() << "); v is " << (Vector&)ball.getVelocity() << "; w is " << (Vector&)ball.getSpin() << "; old state " << ball.getStateString() <<  endl;

    if(ball.getVelocity().mag() > EPSILON
            || u.mag() > EPSILON*ball.getRadius())  //The ball is moving
    {
        if(u.mag() > EPSILON*ball.getRadius())
            newState = Ball::SLIDING;
        else {
            newState = Ball::ROLLING;
        }
    }
    else
    {
        if(fabs(ball.getSpin().z) > EPSILON)
            newState = Ball::SPINNING;
        else
            newState = Ball::STATIONARY;
    }

    if (VERBOSE && ball.getState() != newState) {
        Ball dummy;
        dummy.setState(newState);
        cerr << "Ball " << ball.getID() << " transitioned from " << ball.getStateString() << " to " << dummy.getStateString() << endl;
    }
    ball.setState(newState);
}

void Shot::handleBallCollision(Ball &ball1, Ball &ball2)
{
    if(VERBOSE)
    {
        printf("Handling ball/ball collision\n");

        cerr << "Ball 1: (" << ball1.getVelocity().x << ", "
            << ball1.getVelocity().y << ", "
            << ball1.getVelocity().z << "), ("
            << ball1.getSpin().x << ", "
            << ball1.getSpin().y << ", "
            << ball1.getSpin().z << ")" << endl;
        cerr << "Ball 2: (" << ball2.getVelocity().x << ", "
            << ball2.getVelocity().y << ", "
            << ball2.getVelocity().z << "), ("
            << ball2.getSpin().x << ", "
            << ball2.getSpin().y << ", "
            << ball2.getSpin().z << ")" << endl;
    }

    //First, align the two balls to be on the x-axis.
    double phi = Utils::angle(ball1.getPos(), ball2.getPos(), Utils::START_ZERO, Utils::DEGREES);
    ball1.setVelocity(ball1.getVelocity().rotateDeg(-phi));
    ball2.setVelocity(ball2.getVelocity().rotateDeg(-phi));
    ball1.setSpin(ball1.getSpin().rotateDeg(-phi));
    ball2.setSpin(ball2.getSpin().rotateDeg(-phi));


    double mu_B = 0.01;
    double v0x = ball1.getVelocity().x;// - ball2.getVelocity().x;

    Vector vBB; //relative velocity of balls at contact point
    vBB = (ball1.getVelocity() + (ball2.getVelocity()*(-1)) + (Vector(1, 0, 0).cross(ball1.getSpin() + ball2.getSpin())*-0.028575) );
    double vBB_mag = vBB.mag();

    Vector new_obj_v; //referring to ball2
    new_obj_v.y = ball2.getVelocity().y + mu_B*v0x*vBB.y/vBB_mag;
    new_obj_v.x = ball1.getVelocity().x;
    new_obj_v.z = 0;

    Vector new_cue_v; //ball1
    new_cue_v.y = ball1.getVelocity().y - mu_B*v0x*vBB.y/vBB_mag;
    new_cue_v.x = ball2.getVelocity().x;
    new_cue_v.z = 0;

    Vector new_obj_w, new_cue_w;
    //TODO: Test behavior of spinning, rolling, sliding/rolling cases
    if(ball2.getState() == Ball::STATIONARY || ball2.getState() == Ball::SPINNING)
    {
        //Follow Marlow...
        Vector w_to_add = ((Vector(1, 0, 0).cross(vBB))*(-2.5*mu_B*v0x/vBB_mag/Ball::BALL_RADIUS/1.000575));

        new_obj_w = ( ball2.getSpin() + w_to_add );
        new_cue_w = ( ball1.getSpin() + (w_to_add*(-1)) );
    }
    else
    {
        //Apparently Poolfiz just makes the balls roll...
        new_obj_w.x = -new_obj_v.y/Ball::BALL_RADIUS;
        new_obj_w.y = new_obj_v.x/Ball::BALL_RADIUS;
        new_obj_w.z = ball2.getSpin().z;

        new_cue_w.x = -new_cue_v.y/Ball::BALL_RADIUS;
        new_cue_w.y = new_cue_v.x/Ball::BALL_RADIUS;
        new_cue_w.z = ball1.getSpin().z;
    }

    ball1.setVelocity(new_cue_v);
    ball2.setVelocity(new_obj_v);
    ball1.setSpin(new_cue_w);
    ball2.setSpin(new_obj_w);

    //Finally, unrotate the balls to their original collision angle.
    ball1.setVelocity(ball1.getVelocity().rotateDeg(phi));
    ball2.setVelocity(ball2.getVelocity().rotateDeg(phi));
    ball1.setSpin(ball1.getSpin().rotateDeg(phi));
    ball2.setSpin(ball2.getSpin().rotateDeg(phi));

    //This seems to be necessary to match Poolfiz in corner cases, e.g. break shots...
    if(ball1.getVelocity().mag() < 1e-4 && ball1.getSpin().mag() < 1e-3)
    {
        ball1.setVelocity(Vector());
        ball1.setSpin(Vector(0, 0, ball1.getSpin().z));
    }
    if(ball2.getVelocity().mag() < 1e-4 && ball2.getSpin().mag() < 1e-3)
    {
        ball2.setVelocity(Vector());
        ball2.setSpin(Vector(0, 0, ball2.getSpin().z));
    }

    if(VERBOSE)
    {
        cerr << "Ball 1: (" << ball1.getVelocity().x << ", "
            << ball1.getVelocity().y << ", "
            << ball1.getVelocity().z << "), ("
            << ball1.getSpin().x << ", "
            << ball1.getSpin().y << ", "
            << ball1.getSpin().z << ")" << endl;
        cerr << "Ball 2: (" << ball2.getVelocity().x << ", "
            << ball2.getVelocity().y << ", "
            << ball2.getVelocity().z << "), ("
            << ball2.getSpin().x << ", "
            << ball2.getSpin().y << ", "
            << ball2.getSpin().z << ")" << endl;
    }

    handleStateChange(ball1);
    handleStateChange(ball2);
}

//Given a ball adjacent to a rail, updates the ball's velocity and spin from
//before the collision to after it.
void Shot::handleRailCollision(Ball &ball, Table::BoundaryId rail)
{
    if(VERBOSE)
        printf("Handling rail collision\n");

    Vector v = ball.getVelocity();
    Vector newVelocity;

    if(rail == Table::N_RAIL || rail == Table::S_RAIL) //Horizontal rails
        newVelocity = Vector(0.9*v.x, -0.6*v.y, 0);
    else
        newVelocity = Vector(-0.6*v.x, 0.9*v.y, 0);
    ball.setVelocity(newVelocity);

    //Spin is just divided by 10
    ball.setSpin(ball.getSpin()*0.1);

    //The ball will usually be sliding at this point, whatever it was doing before
    handleStateChange(ball);
}

//Pockets a ball that has reached the pocket.
void Shot::handlePocketed(Ball &ball, Table::Pocket pocket)
{
    ball.setState(Table::stateFromPocket(pocket));
    if(VERBOSE)
        cerr << "Ball " << ball.getIDString() << " pocketed in "
            << Table::pocketName(pocket) << endl;

    ball.setVelocity(Vector());
    ball.setSpin(Vector());
}

//Translates cue strike characteristics into changes in the cue ball's velocity
//and spin.
void Shot::handleCueStrike(Ball &ball, const ShotParams &sp)
{
    //if(ball.getState() != STATIONARY)
    //return;  
    double m = BALL_MASS;
    double M = CUE_MASS;
    double cos = gsl_sf_cos(sp.theta * M_PI/180);
    double sin = gsl_sf_sin(sp.theta * M_PI/180);
    double c = fabs(sqrt(Ball::BALL_RADIUS*Ball::BALL_RADIUS - sp.a*sp.a - sp.b*sp.b));
    //TODO: Check that c is sane

    //correctional constant to match poolfiz
    double c1 = 0.999425;

    //F is the force of the cue stick on the ball
    //Equation (1) in Leckie and Greenspan
    double F = (2*m*sp.v) / (1 + (m/M) + c1*(5.0/(2.0*Ball::BALL_RADIUS*Ball::BALL_RADIUS))*(sp.a*sp.a + sp.b*sp.b*cos*cos + c*c*sin*sin - 2*sp.b*c*cos*sin));
    if(VERBOSE)
        printf("F = %g\n", F);
    //Now, set the velocity of the cue ball
    //In i,j,k coordinate system:
    //velocity is (0, -F*cos/m, 0), making no-jumps assumption (2)
    //angular velocity is (1/I)*(-c*F*sin + b*F*cos, a*F*sin, -a*F*cos) (3)

    //The following equations ((2) and (3) from Leckie/Greenspan) are in a
    //ball-centric coordinate system (ijk). j is parallel to the projection of
    //the cue onto the table; k is the same as z; and the ball's center is at
    //(0, 0, 0). v_rotate returns this orientation to that of the table.
    ball.setVelocity(Vector(0, -F*cos/m, 0).rotateDeg(sp.phi + 90.0));
    ball.setSpin((Vector(-c*F*sin + sp.b*F*cos, sp.a*F*sin, -sp.a*F*cos)*(1/I)).rotateDeg(sp.phi + 90.0));

    //Ball is most likely sliding now
    handleStateChange(ball);
}
*/

//Removes all events involving any of the same balls as the given event.
void Shot::removeRelatedEvents(list<Event*> &events, Event *lastEvent)
{
    for(list<Event*>::iterator itr = events.begin(); itr != events.end(); )
    {
        if(!lastEvent || lastEvent->relatedTo(**itr))
        {
            if (*itr == lastEvent) throw "Cannot remove event.";
            delete (*itr);
            itr = events.erase(itr);
        }
        else
            ++itr;
    }
}

/*
   bool Shot::areRelated(Event &e1, Event &e2)
   {
//The default value for events without ball 2 set is 0 == CUE, so be careful
//about whether or not they actually have two balls involved
if (e1.getBall1() == e2.getBall1()) return true;
if (e1.getType() == Event::BALL_COLLISION) {
BallCollisionEvent* ballEvent1 = dynamic_cast<BallCollisionEvent*>(&e1);
if (ballEvent1->getBall2() == e2.getBall1()) return true;
if (e2.getType() == Event::BALL_COLLISION) {
BallCollisionEvent* ballEvent2 = dynamic_cast<BallCollisionEvent*>(&e2);
if (ballEvent1->getBall2() == ballEvent2->getBall2()) return true;
}
} else if (e2.getType() != Event::BALL_COLLISION) {
BallCollisionEvent* ballEvent2 = dynamic_cast<BallCollisionEvent*>(&e2);
if (e1.getBall1() == ballEvent2->getBall2()) return true;
} //else
return false;

return (e1.getBall1() == e2.getBall1()
|| (e1.getType() == Event::BALL_COLLISION && e1.getBall2() == e2.getBall1())
|| (e2.getType() == Event::BALL_COLLISION && (e1.getBall1() == e2.getBall2()
|| (e1.getType() == Event::BALL_COLLISION
&& e1.getBall2() == e2.getBall2()))));

}*/

/*
   void Shot::addFutureEvents(TableState &state, double curTime, list<Event*> &futureEvents)
   {
//Add expected motion-transition events for all moving balls
for(int i = 0; i < state.getNumBalls(); i++)
{
Event *e = nextTransitionEvent(state.getTable(), state.getBall(i), curTime);
if(e != NULL)
futureEvents.push_back(e);
}
//Add expected ball-ball collision events for all pairs of balls
for(int i = 0; i < state.getNumBalls(); i++)
{
if(state.getBall(i).getState() == Ball::SLIDING || state.getBall(i).getState() == Ball::ROLLING)
{
for(int j = 0; j < state.getNumBalls(); j++)
{
if(i != j)
{
Event *e = nextCollisionEvent(state.getTable(), state.getBall(i), state.getBall(j), curTime);
if(e != NULL)
futureEvents.push_back(e);
}
}
}
}
//Add expected ball-rail and ball-pocket events for all moving balls
for(int i = 0; i < state.getNumBalls(); i++)
{
addBoundaryEvents(state.getTable(), state.getBall(i), curTime, futureEvents);
}
}
*/

void Shot::addRelatedFutureEvents(TableState &state, double curTime, list<Event*> &futureEvents, Ball::Type b1, Ball::Type b2/*set<Ball::Type> &balls*/)
{
    //Optimization opportunity: get the "counts" below pre-computed
    //Add expected motion-transition events for moving balls
    //for(set<Ball::Type>::iterator ballItr = balls.begin(); ballItr != balls.end(); ++ballItr)
    //{
    if (b1 != Ball::UNKNOWN_ID) {
        Event *e = nextTransitionEvent(state.getTable(), state.getBall(b1), curTime);
        if(e != NULL)
            futureEvents.push_back(e);
        if(b2 != Ball::UNKNOWN_ID)
        {
            Event *e = nextTransitionEvent(state.getTable(), state.getBall(b2), curTime);
            if(e != NULL)
                futureEvents.push_back(e);
        }
    }
    //}

    //Add expected ball-ball collision events for all pairs of balls
    for(vector<Ball>::iterator i = state.getBegin(); i != state.getEnd(); ++i)
    {
        if (b1 == Ball::UNKNOWN_ID) {
            Event *e = nextTransitionEvent(state.getTable(), *i, curTime);
            if(e != NULL)
                futureEvents.push_back(e);
            addBoundaryEvents(state.getTable(), *i, curTime, futureEvents);
        }
        if(i->getState() == Ball::SLIDING || i->getState() == Ball::ROLLING)
        {
            for(vector<Ball>::iterator j = state.getBegin(); j != state.getEnd(); ++j)
            {
                if(i != j && (b1 == Ball::UNKNOWN_ID || i->getID() == b1 || j->getID() == b1 || i->getID() == b2 || j->getID() == b2)/*(balls.count((Ball::Type)i) > 0 || balls.count((Ball::Type)j) > 0)*/)
                {
                    Event *e = nextCollisionEvent(state.getTable(), *i, *j, curTime);
                    if (e && (((b1 != i->getID() || b2 !=j->getID()) && (b1 !=j->getID() || b2 !=i->getID())) || e->getTime()-curTime > NEAR_FUTURE_EPSILON)) {
                        futureEvents.push_back(e);
                    }
                }
            }
        }
    }
    //Add expected ball-rail and ball-pocket events for moving balls
    //for(set<Ball::Type>::iterator ballItr = balls.begin(); ballItr != balls.end(); ++ballItr)
    //{
    if (b1 != Ball::UNKNOWN_ID) 
        addBoundaryEvents(state.getTable(), state.getBall(b1), curTime, futureEvents);

    if(b2 != Ball::UNKNOWN_ID)
        addBoundaryEvents(state.getTable(), state.getBall(b2), curTime, futureEvents);

    //}
}

Event* Shot::nextTransitionEvent(const Table &table, Ball &ball, double curTime)
{
    if(ball.getState() == Ball::SLIDING)
    {
        //Equation (11) in Leckie/Greenspan
        double u_o = (ball.getVelocity() + (Vector(0, 0, 1).cross(ball.getSpin())*ball.getRadius())).mag();
        double duration = 2*u_o/(7*table.getMuSliding()*Table::g);
        StateChangeEvent* s = new StateChangeEvent(curTime + duration, ball.getID());
        if (VERBOSE) {
            Vector v1 = Vector(0,0,1).cross(ball.getSpin());
            cerr << "(0,0,1) cross spin:" << v1 << endl;
            v1 = v1 * ball.getRadius();
            cerr << "times radius:" << v1 << endl;
            v1 = v1 + ball.getVelocity();
            cerr << "add velocity:" << v1 << endl;
            double d1 = v1.mag();
            cerr << "magnitude is:" << d1 << " or " << u_o << endl;
            cerr << "next:" <<  *s << ball << endl;
        }
        return s;
    }
    else if(ball.getState() == Ball::ROLLING)
    {
        //Equation (15) in Leckie/Greenspan
        double duration = ((ball.getVelocity()).mag())/(table.getMuRolling()*Table::g);
        return new StateChangeEvent(curTime + duration, ball.getID());
    }
    else if(ball.getState() == Ball::SPINNING)
    {
        //Derived from equation (10) in Leckie/Greenspan
        double duration = fabs((2*ball.getRadius()*ball.getSpin().z)/(5*table.getMuSpinning()*Table::g));
        return new StateChangeEvent(curTime + duration, ball.getID());
    }
    else
    {
        return NULL;
    }
}

Event* Shot::nextCollisionEvent(const Table &table, Ball &ball1, Ball &ball2, double curTime)
{ 

    //ball1.updateState(VERBOSE);
    //ball2.updateState(VERBOSE);

    //TODO: Rewrite section below to reuse code, comments

    if(ball2.getState() == Ball::SLIDING || ball2.getState() == Ball::ROLLING)
    {
        //When ball1 has higher type than ball2 and both are moving, ignore, as
        //we've already made this event
        if(ball1.getID() > ball2.getID())
            return NULL;


        //Both balls are moving
        Vector v1 = ball1.getVelocity();
        Vector v2 = ball2.getVelocity();
        Point r1 = ball1.getPos();
        Point r2 = ball2.getPos();

        //cerr << "V1 " << ball1.v.x << "," << ball1.v.y << "," << ball1.v.z << endl;
        //cerr << "V2 " << ball2.v.x << "," << ball2.v.y << "," << ball2.v.z << endl;
        //cerr << "R1 " << ball1.r.x << "," << ball1.r.y << endl;
        //cerr << "R2 " << ball2.r.x << "," << ball2.r.y << endl;

        double mu1 = (ball1.getState() == Ball::SLIDING ? table.getMuSliding() : table.getMuRolling());
        double mu2 = (ball2.getState() == Ball::SLIDING ? table.getMuSliding() : table.getMuRolling());

        //cerr << "mu1=" << mu1 << " mu2=" << mu2 << endl;

        //cerr << "R=" << R << " w2=" << ball2.w << endl;

        //cerr << "ball states " << ball1.state << " " << ball2.getState() << endl;

        //Vector dbgTemp=v_mult(R, v_cross(v_from(0, 0, 1), ball2.w));

        //cerr << "dbgtmp " << dbgTemp << endl;

        Vector u_hat1, u_hat2;
        if(ball1.getState() == Ball::SLIDING)
            u_hat1 = (v1 + (Vector(0, 0, 1).cross(ball1.getSpin())*ball1.getRadius()));
        else
            u_hat1 = v1;

        if(ball2.getState() == Ball::SLIDING)
            u_hat2 = (v2 + (Vector(0, 0, 1).cross(ball2.getSpin())*ball2.getRadius()));
        else
            u_hat2 = v2;

        if (Utils::fequal(u_hat1.mag(), 0) && !Utils::fequal(v1.mag(), 0)) {
            u_hat1=v1;
        }
        u_hat1=u_hat1.norm();
        if (Utils::fequal(u_hat2.mag(), 0) && !Utils::fequal(v2.mag(), 0)) {
            u_hat2=v2;
        }
        u_hat2=u_hat2.norm();

        //cerr << "Uhat1 " << u_hat1.x << "," << u_hat1.y << "," << u_hat1.z << endl;
        //cerr << "Uhat2 " << u_hat2.x << "," << u_hat2.y << "," << u_hat2.z << endl;

        double ax = 0.5*mu1*Table::g*u_hat1.x - 0.5*mu2*Table::g*u_hat2.x;
        double bx = v2.x - v1.x;
        double cx = r2.x - r1.x;
        double ay = 0.5*mu1*Table::g*u_hat1.y - 0.5*mu2*Table::g*u_hat2.y;
        double by = v2.y - v1.y;
        double cy = r2.y - r1.y;

        //cerr << "ax=" << ax << "; bx=" << bx << "; cx=" << cx <<endl;
        //cerr << "ay=" << ay << "; by=" << by << "; cy=" << cy <<endl;

        //vector<gsl_complex> roots;
        double roots[8];
        int status=solveQuartic(roots, cx*cx + cy*cy - 4*ball1.getRadius()*ball1.getRadius(), 2*bx*cx + 2*by*cy,
                bx*bx + 2*ax*cx + by*by + 2*ay*cy,
                2*ax*bx + 2*ay*by, ax*ax + ay*ay);

        double root = (status?0:leastPositiveRealRoot(roots,(ball1.getPos().to_v()-ball2.getPos().to_v()).dot(ball1.getVelocity()-ball2.getVelocity())>0?NEAR_FUTURE_EPSILON:0));

        if (VERBOSE) { 
            cerr << "looking at balls " << ball1.getIDString() << " and " << ball2.getIDString() << endl;
            for (int i=0; i<8; i++) {
                cerr << roots[i] << " ";
            } 
            cerr << endl;
            cerr << "nextCollision root = " << root << endl;
        }

        /*if(root > 0 && root < NEAR_FUTURE_EPSILON) {

        //check if time is in the near future.  if so, make sure the balls are moving towards each other.
        Point old_dr = ball1.getPos() - ball2.getPos();
        Point new_dr = ball1.getPos() + (ball1.getVelocity()*root).to_p() 
        - ball2.getPos() - (ball2.getVelocity()*root).to_p();
        if (old_dr.x*old_dr.x + old_dr.y*old_dr.y > new_dr.x*new_dr.x + new_dr.y+new_dr.y)
        return new BallCollisionEvent(root + curTime, ball1.getID(), ball2.getID());     
        else
        return NULL;
        } else */
        if (root > 0)
            return new BallCollisionEvent(root + curTime, ball1.getID(), ball2.getID());     
        else
            return NULL;
    }
    else if(ball2.getState() == Ball::STATIONARY || ball2.getState() == Ball::SPINNING)
    {
        //ball2 is stationary
        //TODO: Standardize with above
        Vector v1 = ball1.getVelocity();
        Point r1 = ball1.getPos();
        Point r2 = ball2.getPos();
        double mu = (ball1.getState() == Ball::SLIDING ? table.getMuSliding() : table.getMuRolling());
        //u_hat: should be one thing when sliding, another when rolling
        Vector u_hat;
        if(ball1.getState() == Ball::SLIDING)
        {
            //set u_hat to u
            u_hat = (v1 + (Vector(0, 0, 1).cross(ball1.getSpin())*ball1.getRadius()) );
        }
        else
        {
            //set u_hat to v
            u_hat = v1;
        }
        if (Utils::fequal(u_hat.mag(), 0)) {
            if (Utils::fequal(v1.mag(), 0)) return NULL;
            u_hat=v1;
        }
        u_hat = u_hat.norm();

        //We are preparing the coefficients to go into the equation
        //(axt^2+bxt+cx)^2 + (ayt^2+byt+cy)^2 = (2R)^2
        double ax = -0.5*mu*Table::g*u_hat.x;
        double bx = v1.x;
        double cx = r1.x - r2.x;
        double ay = -0.5*mu*Table::g*u_hat.y;
        double by = v1.y;
        double cy = r1.y - r2.y;

        //t^4: ax^2 + ay^2
        //t^3: 2axbx + 2ayby
        //t^2: bx^2 + 2axcx + by^2 + 2aycy
        //t: 2bxcx + 2bycy
        //1: cx^2 + cy^2 - 4R^2
        //vector<gsl_complex> roots;
        double roots[8];
        int status=solveQuartic(roots, cx*cx + cy*cy - 4*ball1.getRadius()*ball1.getRadius(), 2*bx*cx + 2*by*cy,
                bx*bx + 2*ax*cx + by*by + 2*ay*cy,
                2*ax*bx + 2*ay*by, ax*ax + ay*ay);

        double root = (status?0:leastPositiveRealRoot(roots,(ball1.getPos().to_v()-ball2.getPos().to_v()).dot(ball1.getVelocity()-ball2.getVelocity())>0?NEAR_FUTURE_EPSILON:0));

        if (VERBOSE) { 
            cerr << "looking at balls " << ball1.getIDString() << " and " << ball2.getIDString() << endl;
            for (int i=0; i<8; i++) {
                cerr << roots[i] << " ";
            } 
            cerr << endl;
            cerr << "nextCollision root " << root << endl;
        }

        /*
           if(root > 0 && root < NEAR_FUTURE_EPSILON) {
        //check if time is in the near future.  if so, make sure the balls are moving towards each other.
        Point old_dr = ball1.getPos() - ball2.getPos();
        Point new_dr = ball1.getPos() + (ball1.getVelocity()*root).to_p() 
        - ball2.getPos();
        cerr << "old_dr: " << (old_dr.x*old_dr.x + old_dr.y*old_dr.y) 
        << "new_dr: " << (new_dr.x*new_dr.x + new_dr.y+new_dr.y)<< endl;

        if ( > new_dr.x*new_dr.x + new_dr.y+new_dr.y)
        //return new BallCollisionEvent(root + curTime, ball1.getID(), ball2.getID());     
        else
        //return NULL;
        } */
        if(root > 0)
            return new BallCollisionEvent(root + curTime, ball1.getID(), ball2.getID());
        else
            return NULL;
    }
    else
    {
        //Can't hit a ball that's not on the table
        return NULL;
    }
}

void Shot::addBoundaryEvents(const Table &table, Ball &ball, double curTime, list<Event*> &futureEvents)
{
    //TODO:check rail collisions (at least), some balls are going off the table when I test them.
    if(!(ball.getState() == Ball::SLIDING || ball.getState() == Ball::ROLLING))
        return;

    //Rail boundary events
    int numRoots;
    double root1, root2; //Store the results of solving the quadratic
    double eventTime;
    double collision_y;
    double mu;
    double t;
    if(ball.getState() == Ball::SLIDING)
        mu = table.getMuSliding();
    else
        mu = table.getMuRolling();
    Point r = ball.getPos();
    Vector v = ball.getVelocity();
    Vector u_hat;
    if(ball.getState() == Ball::SLIDING)
        u_hat = (v + (Vector(0, 0, 1).cross(ball.getSpin())*ball.getRadius())).norm();
    else
        u_hat = v.norm(); //When rolling, use v-hat instead
    double cp_width = table.getPocketLeft(Table::SW).x; //width of a corner pocket over sqrt(2)
    double t_length = table.getLength();
    double t_width = table.getWidth();

    //Equations based on suggestions of Leckie/Greenspan, section 3.2.2.
    //TODO: Consider decomping this section, to make it less ugly

    //South rail: y=0
    numRoots = gsl_poly_solve_quadratic(-0.5*mu*Table::g*u_hat.y, v.y, r.y - ball.getRadius(), &root1, &root2);
    eventTime = calcEventTime(numRoots, root1, root2, curTime);
    if(Utils::fgreater(eventTime, curTime)) //eventTime > curTime + EPSILON)
        futureEvents.push_back(new RailCollisionEvent(eventTime, ball.getID(), Table::S_RAIL));

    //North rail: y=t_length
    numRoots = gsl_poly_solve_quadratic(-0.5*mu*Table::g*u_hat.y, v.y, r.y - t_length + ball.getRadius(), &root1, &root2);
    eventTime = calcEventTime(numRoots, root1, root2, curTime);
    if(Utils::fgreater(eventTime, curTime)) //eventTime > curTime + EPSILON)
        futureEvents.push_back(new RailCollisionEvent(eventTime, ball.getID(), Table::N_RAIL));

    //West rails and side pocket: x=0
    numRoots = gsl_poly_solve_quadratic(-0.5*mu*Table::g*u_hat.x, v.x, r.x - ball.getRadius(), &root1, &root2);
    eventTime = calcEventTime(numRoots, root1, root2, curTime);
    if(Utils::fgreater(eventTime, curTime)) //eventTime > curTime + EPSILON)
    {
        t = eventTime - curTime;
        collision_y = r.y + v.y*t - 0.5*mu*Table::g*u_hat.y*t*t;
        if(collision_y < table.getPocketLeft(Table::W).y)
            futureEvents.push_back(new RailCollisionEvent(eventTime, ball.getID(), Table::SW_RAIL));
        else if(collision_y > table.getPocketRight(Table::W).y)
            futureEvents.push_back(new RailCollisionEvent(eventTime, ball.getID(), Table::NW_RAIL));
        else
            futureEvents.push_back(new PocketedEvent(eventTime, ball.getID(), Table::W));
    }

    //East rails and side pocket: x=t_width
    numRoots = gsl_poly_solve_quadratic(-0.5*mu*Table::g*u_hat.x, v.x, r.x - t_width + ball.getRadius(), &root1, &root2);
    eventTime = calcEventTime(numRoots, root1, root2, curTime);
    if(Utils::fgreater(eventTime, curTime)) //eventTime > curTime + EPSILON)
    {
        t = eventTime - curTime;
        collision_y = r.y + v.y*t - 0.5*mu*Table::g*u_hat.y*t*t;
        if(collision_y < table.getPocketRight(Table::E).y)
            futureEvents.push_back(new RailCollisionEvent(eventTime, ball.getID(), Table::SE_RAIL));
        else if(collision_y > table.getPocketLeft(Table::E).y)
            futureEvents.push_back(new RailCollisionEvent(eventTime, ball.getID(), Table::NE_RAIL));
        else
            futureEvents.push_back(new PocketedEvent(eventTime, ball.getID(), Table::E));
    }

    //SW pocket
    numRoots = gsl_poly_solve_quadratic(-0.5*mu*Table::g*(u_hat.y + u_hat.x), v.x + v.y, r.x + r.y - cp_width - ball.getRadius(), &root1, &root2);
    eventTime = calcEventTime(numRoots, root1, root2, curTime);
    if(eventTime >= 0)
        futureEvents.push_back(new PocketedEvent(eventTime, ball.getID(), Table::SW));

    //SE pocket
    numRoots = gsl_poly_solve_quadratic(0.5*mu*Table::g*(u_hat.y - u_hat.x), v.x - v.y, r.x - r.y + cp_width + ball.getRadius() - t_width, &root1, &root2);
    eventTime = calcEventTime(numRoots, root1, root2, curTime);
    if(eventTime >= 0)
        futureEvents.push_back(new PocketedEvent(eventTime, ball.getID(), Table::SE));

    //NW pocket
    numRoots = gsl_poly_solve_quadratic(0.5*mu*Table::g*(u_hat.x - u_hat.y), v.y - v.x, r.y - r.x + cp_width + ball.getRadius() - t_length, &root1, &root2);
    eventTime = calcEventTime(numRoots, root1, root2, curTime);
    if(eventTime >= 0)
        futureEvents.push_back(new PocketedEvent(eventTime, ball.getID(), Table::NW));

    //NE pocket
    numRoots = gsl_poly_solve_quadratic(-0.5*mu*Table::g*(u_hat.x + u_hat.y), v.x + v.y, r.x + r.y + cp_width + ball.getRadius() - t_width - t_length, &root1, &root2);
    eventTime = calcEventTime(numRoots, root1, root2, curTime);
    if(eventTime >= 0)
        futureEvents.push_back(new PocketedEvent(eventTime, ball.getID(), Table::NE));
}

//Returns the smaller non-negative root plus curTime. If there are no non-negative
//roots, returns -1.
double Shot::calcEventTime(int numRoots, double root1, double root2, double curTime)
{
    switch(numRoots)
    {
        case 0:
            return -1;
        case 1:
            return (root1 >= 0 ? root1 + curTime : -1);
        case 2:
            if(root1 >= 0)
            {
                if(root2 >= 0)
                    return (root1 < root2 ? root1 + curTime : root2 + curTime);
                else
                    return root1 + curTime;
            }
            else
            {
                if(root2 >= 0)
                    return root2 + curTime;
                else
                    return -1;
            }
        default:
            printf("Error: wrong input %d to Shot::calcEventTime\n", numRoots);
            return -2;
    }
}

int Shot::solveQuartic(double roots[], double a0, double a1, double a2, double a3, double a4)
{
    double input[5];
    input[0] = a0;
    input[1] = a1;
    input[2] = a2;
    input[3] = a3;
    input[4] = a4;

    gsl_poly_complex_workspace *w = gsl_poly_complex_workspace_alloc(5);
    if(w == NULL)
    {
        printf("Error: gsl workspace not allocated in Shot::solveQuartic\n");
        return 1;
    }
    int status=gsl_poly_complex_solve (input, 5, w, roots);
    /*if(gsl_poly_complex_solve (input, 5, w, results) != GSL_SUCCESS)
      {
    //Probably not really a problem...
    printf("Warning: quartic wasn't solved in Shot::solveQuartic\n");
    return;
    }*/

    gsl_poly_complex_workspace_free(w);
    return status;
}

double Shot::leastPositiveRealRoot(double roots[], double epsilon)
{
    double leastRoot = -1;

    for(int i = 0; i < 4; i++)
    {
        if(Utils::fequal(roots[i*2+1], 0) && roots[i*2]> epsilon) //Ignore complex or negative roots, or 0
        {
            if(leastRoot < 0 || roots[i*2] < leastRoot)
                leastRoot = roots[i*2];
        }
    }

    return leastRoot;
}

ostream& operator<<(ostream &out, TableState &rhs) {
    out << "TableState {";
    for (vector<Ball>::iterator i=rhs.getBegin(); i != rhs.getEnd(); ++i) {
        out << endl << " " << *i;
    }
    out << endl << "}";
    return out;
}

TableState& TableState::operator=(const TableState &rhs) {
    if (this != &rhs) {
        if (&_table != &rhs._table) throw "Incompatible TableState assignment.";
        balls = rhs.balls;
    }
    return *this;
}

void TableState::setBall(Ball::Type btype, Ball::State state, Point r) {
    if ((int)btype > (int)Ball::UNKNOWN_ID) 
        throw "Illegal Ball::Type to get from TableState"; // + (int)btype; 
    for (vector<Ball>::iterator iter = balls.begin(); iter != balls.end(); ++iter) {
        if (iter->getID() == btype) {
            iter->setState(state);
            iter->setPos(r);
            return;
        }
    } //else
    balls.push_back(Ball(btype, state, r));
}

void TableState::setBall(Ball &b) {
    setBall(b.getID(), b.getState(), b.getPos());
}

void TableState::setBall(Ball::Type btype, Ball::State state, double x, double y){
    setBall(btype, state, Point(x, y));
} 

void TableState::spotBall(Ball::Type b, double dither) {
    setBall(b, Ball::STATIONARY, _table.getFootSpot());
    getBall(b).addNoise(dither);
    int other = findOverlap(b);
    int tries = 0;
    //cerr << "0th try: " << other << " " << getBall(b) << endl;
    while (other > -1 && tries < 100) {      
        tries++;

        double dx = balls[other].getPos().x - _table.getWidth()/2;
        double dr = balls[other].getRadius() + getBall(b).getRadius();
        double newy = balls[other].getPos().y - sqrt(dr*dr - dx*dx);
        setBall(b, Ball::STATIONARY, _table.getWidth()/2, 
                newy);
        getBall(b).addNoise(dither);
        newy = newy - fabs(newy - getBall(b).getPos().y);
        if (newy < getBall(b).getRadius())
            throw "Could not spot ball, ball out of bounds";
        setBall(b, Ball::STATIONARY, getBall(b).getPos().x, 
                newy);
        //cerr << "diff is" << (balls[other].getPos().y - 
        //                      getBall(b).getPos().y - 
        //                      balls[other].getRadius() - 
        //                      getBall(b).getRadius()) << endl;
        other = findOverlap(b);
        //cerr << tries << "th try: " << other << " " << getBall(b) << endl;
    }
    if (other > -1) {
        throw "Could not spot ball"; // + getBall(b).getIDString();
    }
}

int TableState::findOverlap(Ball::Type ball) const {
    for (unsigned int i=0; i<balls.size(); i++) {
        if (balls[i].getID() != ball && balls[i].overlaps(getBall(ball))) 
            return i;
    }
    return -1;
}

const Ball& TableState::getBall(Ball::Type btype) const {
    return const_cast<const Ball&>(const_cast<TableState*>(this)->getBall(btype));
}

Ball& TableState::getBall(Ball::Type btype) {
    if ((int)btype > (int)Ball::UNKNOWN_ID) 
        throw "Illegal Ball::Type to get from TableState"; // + (int)btype; 

    vector<Ball>::iterator iter;
    for (iter = balls.begin(); iter != balls.end(); ++iter)
        if (iter->getID() == btype) return *iter;
    balls.push_back(Ball(btype, Ball::NOTINPLAY));
    return balls.back();
}

int TableState::isValidBallPlacement(bool VERBOSE) const
{
    int result = OK_PRECONDITION;
    if (getNumBalls() < 1) return result; //There are no balls, so it's automatically valid

    /*if (ignore2<ignore1) {
      Ball::Type tmp=ignore1;
      ignore1=ignore2;
      ignore2=tmp;
      }*/

    // look for any ball off the table
    double ball1r = balls[0].getRadius();
    for (vector<Ball>::const_iterator ball1 = balls.begin(); ball1 != balls.end(); ++ball1) {
        //ball1r = s.getBall(ball1).getRadius(); Uncomment if balls are expected to have different radii
        if (ball1->isInPlay()) {
            if ( Utils::fless(ball1->getPos().x, ball1r) || Utils::fgreater(ball1->getPos().x, _table.getWidth()-ball1r) ) {
                result |= BAD_X_VAL;
                if (VERBOSE) cerr << "bad x val: " << *ball1 << endl;
            }
            if ( Utils::fless(ball1->getPos().y, ball1r) || Utils::fgreater(ball1->getPos().y, _table.getLength()-ball1r) ) { 
                result |= BAD_Y_VAL;
                if (VERBOSE) cerr << "bad y val: " << *ball1 << endl;
            }
        }
    }

    // look for any ball overlap with all other balls
    //double ball2r = ball1r;
    bool overlap = false;
    for (vector<Ball>::const_iterator ball1 = balls.begin(); ball1 != balls.end()-1; ++ball1) {
        if (overlap) break;
        if (ball1->isInPlay()) {
            //ball1r = s.ball1->getRadius(); Uncomment if balls are expected to have different radii
            for (vector<Ball>::const_iterator ball2 = ball1+1; ball2 != balls.end(); ++ball2) {
                if (ball2->isInPlay()) {
                    //ball2r = s.ball2->getRadius(); Uncomment if balls are expected to have different radii
                    if (ball1->overlaps(*ball2) /* && (ball1->getID()!=ignore1 || ball2->getID()!=ignore2)*/) {
                        if (VERBOSE) {

                            cerr << "overlap: " << *ball1 << " and " << *ball2 << endl;
                            cerr << setprecision(20) << (ball1->dist(*ball2)/2) << endl;
                            cerr << setprecision(6);
                        }
                        result |= BALL_OVERLAP;
                        overlap = true;
                        break;
                    } 	
                } 
            }
        } 
    }
    return result;
}

int TableState::isPhysicallyPossible(const ShotParams &sp,bool VERBOSE) const
{
    int result = OK_PRECONDITION;
    double a = sp.a;
    double b = sp.b;
    double theta = sp.theta;
    double phi = sp.phi;
    double V = sp.v;

    // add epsilon to b and theta to prevent artificial CUE_STICK_COLLISION due to numerical precision on different machines
    theta += EPSILON_THETA;	
    b += EPSILON_B;

    theta *= M_PI/180;	// convert to radians
    phi *= M_PI/180;	// convert to radians
    a /= 1000;			// convert to metres
    b /= 1000;			// convert to metres
    double a2 = pow(a,2);
    double b2 = pow(b,2);		
    double rad = getBall(Ball::CUE).getRadius();

    // check shot parameters
    double c = rad*rad - a2 - b2;
    if ( Utils::fgreater(V, MAX_VELOCITY) || Utils::fless(V, 0.0) ) result |= BAD_V_VAL;
    if ( Utils::fgreater(fabs(phi), 2*M_PI) ) result |= BAD_PHI_VAL;

    double thetaMax;
    thetaMax = MAX_THETA*M_PI/180.0;
    if ( Utils::fgreater(theta, thetaMax) || Utils::fless(theta, MIN_THETA) ) result |= BAD_THETA_VAL;
    if ( Utils::fless(fabs(a), rad) && Utils::fless(fabs(b), rad) && Utils::fless(c, 0.0) ) {
        result |= BAD_A_VAL;
        result |= BAD_B_VAL;
    }
    if ( Utils::fgreater(fabs(a), rad) ) result |= BAD_A_VAL;
    if ( Utils::fgreater(fabs(b), rad) ) result |= BAD_B_VAL;

    // check for cue-ball collisions:
    c = sqrt(c);

    // coords of impact point in ball frame
    Vector p1(a, c, b);  

    // coords in table frame
    Point rb = getBall(Ball::CUE).getPos();
    Vector r(rb.x, rb.y, rad);
    p1 = r + p1.rotateRad(-1*(phi + M_PI/2));

    // unit vector in table plane from impact point along cue direction
    Vector q( cos(phi), sin(phi), 0 );
    q = q.norm();

    // coords of another point on the projection of the cue in the table frame
    Vector pi = p1 - (q * (_table.getCueLength() * cos(theta)));

    // coords of the butt end of the cue in the table frame
    Vector p2( pi.x, pi.y, pi.z + _table.getCueLength() * sin(theta) );

    //Vector pcue(state.getBall(CUE).getPos().x, state.getBall(CUE).getPos().y, rad);

    // for all balls, check for intersections with cue
    for (vector<Ball>::const_iterator i = balls.begin(); i != balls.end(); ++i) {
        if (i->isInPlay() && i->getID() != Ball::CUE) {
            Vector p3( i->getPos().x, i->getPos().y, rad );	
            //Vector intr1, intr2;
            double root1, root2;
            int numIntr = numLineSphereIntersections(p1, p2, p3, rad, root1, root2);
            if (numIntr > 0 && 
                    ((Utils::fgreater(root1, 0) && Utils::fless(root1,_table.getCueLength())) ||
                     (Utils::fgreater(root2, 0) && Utils::fless(root2,_table.getCueLength()))))
            {
                result |= CUE_STICK_COLLISION;
                if (VERBOSE) {
                    cerr<<"Cue collided with ball " << i->getIDString() << ", where CUE is " << p1 << " to " << p2 << " and ball is " << p3 << " " << rad << endl;
                }
                break;
            }
        }			
    } // end for all balls

    // check for collisions between cue stick and rails:
    double t, h;
    double testHeight = _table.getRailHeight();
    int railIdx = 0;
    while ( ((result & CUE_STICK_COLLISION) == 0) && railIdx < 4 ) {
        switch (railIdx) {
            case 0: 	
                if ( Utils::fgreater(p1.y, 0.0) && Utils::fless(p2.y, 0.0) ) { // cue crosses S rail
                    t = (0.0 - p1.y) / (p2.y - p1.y);
                    h = p1.z + t*(p2.z - p1.z);	
                    if ( Utils::flessequal(h, testHeight) ) {
                        result |= CUE_STICK_COLLISION;			
                        if (VERBOSE) {
                            cerr<<"Cue collided with S rail CUE is " << p1 << p2 << " t = " <<t << " h = " <<h <<endl;
                        }
                    }
                }
                break;
            case 1: 	
                if ( Utils::fless(p1.x, _table.getWidth()) && Utils::fgreater(p2.x, _table.getWidth()) ) { // cue crosses NE or SE rail
                    t = (_table.getWidth() - p1.x) / (p2.x - p1.x);
                    h = p1.z + t*(p2.z - p1.z);	
                    if ( Utils::flessequal(h, testHeight) ) {
                        result |= CUE_STICK_COLLISION;			
                        if (VERBOSE) {
                            cerr<<"Cue collided with E rail CUE is " << p1 << p2 << " t = " <<t << " h = " <<h <<endl;
                        }
                    }
                }
                break;
            case 2: 	
                if ( Utils::fgreater(p1.x, 0.0) && Utils::fless(p2.x, 0.0) ) { // cue crosses NW or SW rail
                    t = (0.0 - p1.x) / (p2.x - p1.x);
                    h = p1.z + t*(p2.z - p1.z);	
                    if ( Utils::flessequal(h, testHeight) ) {
                        result |= CUE_STICK_COLLISION;			
                        if (VERBOSE) {
                            cerr<<"Cue collided with W rail CUE is " << p1 << p2 << " t = " <<t << " h = " <<h <<endl;
                        }
                    }
                }
                break;
            case 3: 	
                if ( Utils::fless(p1.y, _table.getLength()) && Utils::fgreater(p2.y, _table.getLength()) ) { // cue crosses N rail
                    t = (_table.getLength() - p1.y) / (p2.y - p1.y);
                    h = p1.z + t*(p2.z - p1.z);	
                    if ( Utils::flessequal(h, testHeight) ) {
                        result |= CUE_STICK_COLLISION;			
                        if (VERBOSE) {
                            cerr<<"Cue collided with N rail CUE is " << p1 << p2 << " t = " <<t << " h = " <<h <<endl;
                        }
                    }
                }
                break;
        }  
        railIdx++;		
    } // end while

    return result;
}

void TableState::addNoise(double dither)
{
    //no balls to get noisy positions for
    if (balls.empty()) return;

    //balls are not in a valid position
    if (isValidBallPlacement() != OK_PRECONDITION) {
        throw "Balls are not in a valid position before adding noise";
    }
    // set up an array of noisy positions
    TableState noisypostns(_table);
    bool rackDone = false;
    int rackTries = 0;
    while (!rackDone && rackTries < 10000) {
        rackTries++;

        //copy the balls over
        noisypostns.balls = balls;

        // dither the non-noisy positions
        for (vector<Ball>::iterator i = noisypostns.getBegin(); i != noisypostns.getEnd(); ++i) {
            i->addNoise(dither);
        }

        // check for overlap
        rackDone = (noisypostns.isValidBallPlacement() == OK_PRECONDITION);
    }

    // copy the noisy positions back to the original array if noisifying was successful
    if (rackDone) {
        balls = noisypostns.balls;
    } else {
        throw "Could not add noise to balls";
    }
}

void Ball::addNoise(double dither) {
    // set up random number generator
    double rngSigma = radius * dither / 3.0;
    //cerr << "old:" << r;
    r.x += gsl_ran_gaussian(Utils::rng(), rngSigma);
    r.y += gsl_ran_gaussian(Utils::rng(), rngSigma);
    //cerr << "new:" << r << endl;
}

void TableState::randomize() {
    //no balls to randomize
    if (balls.empty()) return;

    TableState randomBalls(*this);

    srand(time(0));
    for (vector<Ball>::iterator b = randomBalls.getBegin(); b != randomBalls.getEnd(); ++b) {
        if (b->getState() != Ball::NOTINPLAY) {
            int tries = 0;
            bool overlap = false;
            do {
                overlap = false;
                tries ++;
                b->setPos(Point(( gsl_rng_uniform(Utils::rng()) * 
                                (_table.getWidth() - (b->getRadius()*2)) ) + b->getRadius(), 
                            ( gsl_rng_uniform(Utils::rng()) * 
                              (_table.getLength() - (b->getRadius()*2)) ) + b->getRadius()));
                for (vector<Ball>::iterator b2 = randomBalls.getBegin(); (b2 < b) && (overlap == false); ++b2) {
                    if (b2->getState() != Ball::NOTINPLAY) {
                        //cerr << (b2 < b) << "b is" << *b << "b2 is" << *b2 << endl;
                        overlap = b->overlaps(*b2);
                    }
                }
                //cerr << "overlap:" << overlap << " tries:" << tries << endl;
            } while (overlap && tries < 100);
            if (overlap) {
                throw "Could not randomize TableState";
            }
            b->setState(Ball::STATIONARY);
        }
    }
    balls = randomBalls.balls;

}

Ball::Type TableState::getFirstBallHit(const ShotParams &sp) {
    return Shot(*this, sp, false).getEventList().back()->getBall2();
}

void TableState::toStream(ostream &out) const {
    out << balls.size() << " ";
    for (unsigned int i=0; i<balls.size(); i++) {
        balls[i].toStream(out);
    }
}

void TableState::fromStream(istream &in) {
    int numBalls;
    in >> numBalls;
    if (numBalls>MAX_BALLS_ON_TABLE) throw "Too many balls on table!\n";
    Ball newBall;
    for (int i=0; i<numBalls; i++) {
        newBall.fromStream(in);
        balls.push_back(newBall);
    }
}

string TableState::toString() const {
    ostringstream out;
    toStream(out);
    return out.str();
}

void TableState::fromString(const string &s) {
    istringstream in(s);
    fromStream(in);
}

int TableState::numLineSphereIntersections(Vector &p1, Vector &p2, Vector &p3, double rad, double& root1, double& root2) {
    root1 = -1;
    root2 = -1;
    Vector l = (p2 - p1).norm();
    Vector c(p3 - p1);
    double det = pow(l.dot(c), 2.0) - c.dot(c) + pow(rad, 2.0);
    if ( Utils::fless(det, 0.0) ) return 0;
    else if ( Utils::fgreater(det, 0.0) ) {
        root1 = l.dot(c)-det;
        root2 = l.dot(c)+det;
        return 2;
    }
    else {
        root1 = l.dot(c);
        return 1;
    }
}

string getFastFizVersion() {
    return "FastFiz version 0.1 built " __DATE__ " " __TIME__;
}

void rack(TableState& ts) {

    /*
       const Table& t = ts.getTable();
       double w = t.getWidth();
       double x = t.getFootSpot().x;
       double y = t.getFootSpot().y;
       double r = Ball::BALL_RADIUS;
       double e = r;

       double dy = e/2*sqrt(3)+r*sqrt(3);
       double dx = e/2+r;


       ts.setBall(Ball::CUE,      Ball::STATIONARY, w/2, t.getHeadString());
       ts.setBall(Ball::ONE,      Ball::STATIONARY, x, y);
    //ts.getBall(Ball::ONE).setState(Ball::SLIDING);
    //ts.getBall(Ball::ONE).setSpin(Vector(180,0,0));

    ts.setBall(Ball::TWO,      Ball::STATIONARY, x-dx,   y-dy); 
    ts.setBall(Ball::THREE,    Ball::STATIONARY, x+dx,   y-dy); 
    ts.setBall(Ball::FOUR,     Ball::STATIONARY, x-dx*2, y-dy*2); 
    ts.setBall(Ball::FIVE,     Ball::STATIONARY, x,      y-dy*2); 
    ts.setBall(Ball::SIX,      Ball::STATIONARY, x+dx*2, y-dy*2); 
    ts.setBall(Ball::SEVEN,    Ball::STATIONARY, x-dx*3, y-dy*3); 
    ts.setBall(Ball::EIGHT,    Ball::STATIONARY, x-dx,   y-dy*3); 
    ts.setBall(Ball::NINE,     Ball::STATIONARY, x+dx,   y-dy*3); 
    ts.setBall(Ball::TEN,      Ball::STATIONARY, x+dx*3, y-dy*3); 
    ts.setBall(Ball::ELEVEN,   Ball::STATIONARY, x-dx*4, y-dy*4); 
    ts.setBall(Ball::TWELVE,   Ball::STATIONARY, x-dx*2, y-dy*4); 
    ts.setBall(Ball::THIRTEEN, Ball::STATIONARY, x,      y-dy*4); 
    ts.setBall(Ball::FOURTEEN, Ball::STATIONARY, x+dx*2, y-dy*4); 
    ts.setBall(Ball::FIFTEEN,  Ball::STATIONARY, x+dx*4, y-dy*4); 

    ts.addNoise(10);
    */

    //ts.fromString("16 0.028575 1 0 0.37744090581937889128 2.1191419404969584761 0.028574999999999999706 1 1 0.55499705417823386178 0.6646260604392125737 0.028574999999999999706 1 2 0.68724784068010624782 0.65243101981344653328 0.028574999999999999706 1 3 0.60219446359108330658 0.44304506094407225536 0.028574999999999999706 1 4 0.41672000254982610734 0.60385126922524734461 0.028574999999999999706 1 5 0.63677871238113648023 0.52010138332172162201 0.028574999999999999706 1 6 0.48862715017378310911 0.11097820143256829917 0.028574999999999999706 1 7 0.31247768202523606984 0.41713155216634056899 0.028574999999999999706 1 8 0.4995072322551756816 0.51718984437937587373 0.028574999999999999706 1 9 0.78060556467038433315 0.40973949569223805378 0.028574999999999999706 1 10 0.87885070403890341861 0.26055420980464616409 0.028574999999999999706 1 11 0.39833791352033265376 0.15362918052031551697 0.028574999999999999706 1 12 0.32007028290170003171 0.27575892620947933809 0.028574999999999999706 1 13 0.42942861921247355017 0.34812774558809372882 0.028574999999999999706 1 14 0.65937202467564970387 0.1249417008847652083 0.028574999999999999706 1 15 0.53391732180295847776 0.2401589277373836484 ");

    //ts.fromString("16 0.028575 1 0 0.67979039293473186856 1.9138650440890030424 0.028574999999999999706 1 1 0.57538246007141813365 0.41287332516915370428 0.028574999999999999706 1 2 0.53885089033550614968 0.66989194634159421327 0.028574999999999999706 1 3 0.52172485394339396247 0.47830230992545896829 0.028574999999999999706 1 4 0.33725948094313989989 0.26437125516220322252 0.028574999999999999706 1 5 0.59596930567740424856 0.49171371369358712888 0.028574999999999999706 1 6 0.67035935540182745029 0.43507334154379884161 0.028574999999999999706 1 7 0.48515307980114197317 0.37044648547624714485 0.028574999999999999706 1 8 0.44990141093053021004 0.26574422065492098177 0.028574999999999999706 1 9 0.55336656098679071203 0.32680582129985524809 0.028574999999999999706 1 10 0.62144754826452308638 0.36269604598132215711 0.028574999999999999706 1 11 0.40167784427735303554 0.21427754466404583256 0.028574999999999999706 1 12 0.63331157847006702699 0.28483379611254233676 0.028574999999999999706 1 13 0.52082679209380655649 0.20673464698058210964 0.028574999999999999706 1 14 0.49384237605007397232 0.081963882980097629849 0.028574999999999999706 1 15 0.61333720327411744844 0.18175695237011210703 ");

    //ts.fromString("16 0.028575 1 0 0.70330654119767865851 1.7321769298305966789 0.028574999999999999706 1 1 0.67705867019634091708 0.78956972363215383037 0.028574999999999999706 1 2 0.63201432786743061421 0.61593901698959674373 0.028574999999999999706 1 3 0.82523721783183690004 0.574192112153682932 0.028574999999999999706 1 4 0.3699620464274247511 0.31745000561861475008 0.028574999999999999706 1 5 0.43321051203452926037 0.37340818068060704915 0.028574999999999999706 1 6 0.66482751949097151911 0.43863119024119451295 0.028574999999999999706 1 7 0.42486065661746558186 0.4386771638717514854 0.028574999999999999706 1 8 0.52181881387708395348 0.45960060339713282307 0.028574999999999999706 1 9 0.87319745530365011632 0.25873797647734364524 0.028574999999999999706 1 10 0.79246755789084410981 0.30493935063838600863 0.028574999999999999706 1 11 0.26231163944132090959 0.18624583265861593961 0.028574999999999999706 1 12 0.5610639762518883078 0.33466667602296368278 0.028574999999999999706 1 13 0.73317637613445740818 0.21695194732136688365 0.028574999999999999706 1 14 0.59072710343288381551 0.28119927934717242612 0.028574999999999999706 1 15 0.76020010854742114947 0.16545217035060477695 ");

    //ts.fromString("16 0.028575 1 0 0.56550967282332587072 1.7311309235144776153 0.028574999999999999706 1 1 0.58957693944102196326 0.91228190489269844754 0.028574999999999999706 1 2 0.43374334935826580617 0.47156595720319655074 0.028574999999999999706 1 3 0.50078425397356773896 0.4734060659668859472 0.028574999999999999706 1 4 0.6195653528980203939 0.52491846925800844659 0.028574999999999999706 1 5 0.47236410635454767482 0.2994309689920531814 0.028574999999999999706 1 6 0.4945407817509851478 0.37605639535101226256 0.028574999999999999706 1 7 0.23876000672331212926 0.41828187095050678446 0.028574999999999999706 1 8 0.67255309314892408601 0.32769005618420149473 0.028574999999999999706 1 9 0.54112747730520449618 0.32366129590608372357 0.028574999999999999706 1 10 0.7276183273735344903 0.26534174272505400838 0.028574999999999999706 1 11 0.50038616782329314869 0.0975085923432981172 0.028574999999999999706 1 12 0.44107853186553047653 0.23368483004931958624 0.028574999999999999706 1 13 0.36148599890668231538 0.14230747893400674808 0.028574999999999999706 1 14 0.62997138459244317321 0.17044871426330535091 0.028574999999999999706 1 15 0.75576637154205739133 0.2125264686846653206 ");

    //ts.fromString("16 0.028575 1 0 0.50379200881412056212 1.7343431930493833359 0.028574999999999999706 1 1 0.58973625913905758367 0.64106884877716752147 0.028574999999999999706 1 2 0.48430438570990780311 0.64295228791120140333 0.028574999999999999706 1 3 0.76832245754271488725 0.52367020811710152106 0.028574999999999999706 1 4 0.43475607224147921892 0.43951428483891136212 0.028574999999999999706 1 5 0.65730255233019607708 0.31855616197511577026 0.028574999999999999706 1 6 0.82561406775136803038 0.31219564038072039613 0.028574999999999999706 1 7 0.46849531916615655414 0.127617293876370691 0.028574999999999999706 1 8 0.52651016693734664642 0.21136621139421071791 0.028574999999999999706 1 9 0.5876921229525807755 0.27170248274445291115 0.028574999999999999706 1 10 0.65945622351337107592 0.38447178305394635878 0.028574999999999999706 1 11 0.3099839963492645456 0.2136203000299724819 0.028574999999999999706 1 12 0.49895880124494929042 0.4529839523175398508 0.028574999999999999706 1 13 0.65002840204137968261 0.18455037376373220614 0.028574999999999999706 1 14 0.69136070329023713299 0.26440597250057906731 0.028574999999999999706 1 15 0.87894055027614326203 0.13360566963197900714");

    //ts.fromString("16 0.028575 1 0 0.68722594436896677905 1.7040906040187280279 0.028574999999999999706 1 1 0.53382050422751237129 0.63162008402187275191 0.028574999999999999706 1 2 0.46772151323514554733 0.61587134576935986807 0.028574999999999999706 1 3 0.72337637967766998948 0.53344628975923358283 0.028574999999999999706 1 4 0.46477350666674693835 0.25267731816514227283 0.028574999999999999706 1 5 0.35100696623023802667 0.4383650664070787295 0.028574999999999999706 1 6 0.66165239900005679008 0.5979681430893351024 0.028574999999999999706 1 7 0.40975087682933203359 0.31440964448100955364 0.028574999999999999706 1 8 0.49866886439917640361 0.39434222197876672711 0.028574999999999999706 1 9 0.60141813631477203383 0.47415994892849228082 0.028574999999999999706 1 10 0.58466870557343264547 0.36154030497665351485 0.028574999999999999706 1 11 0.25447667697032294587 0.22625074628058805537 0.028574999999999999706 1 12 0.49940341852224606756 0.11139314895829204033 0.028574999999999999706 1 13 0.42049422337565944829 0.14611081055048505406 0.028574999999999999706 1 14 0.6979610963871016871 0.14492696313348665815 0.028574999999999999706 1 15 0.6767771654429949324 0.29464052542703045301 ");

    //ts.fromString("16 0.028575 1 0 0.3855646249752190946 1.5588653457407597447 0.028574999999999999706 1 1 0.62326601583864615908 0.68272937617348616879 0.028574999999999999706 1 2 0.43271043577117307422 0.45335376340314509358 0.028574999999999999706 1 3 0.51613308084129305708 0.60365544236064472283 0.028574999999999999706 1 4 0.36139279957327558046 0.57305510951936411423 0.028574999999999999706 1 5 0.67253129297124969455 0.64074505118128444181 0.028574999999999999706 1 6 0.69221884944746325807 0.34443177444473177529 0.028574999999999999706 1 7 0.5128737578684485765 0.26815379698731273139 0.028574999999999999706 1 8 0.53877189338046871381 0.50158783702227915935 0.028574999999999999706 1 9 0.62452591412144919936 0.38331593030323474203 0.028574999999999999706 1 10 0.73919169154631114704 0.46411848838130459471 0.028574999999999999706 1 11 0.42261868544269481918 0.07262065605147470515 0.028574999999999999706 1 12 0.47262425841082783551 0.34323853175928337178 0.028574999999999999706 1 13 0.60027058631873664307 0.30857765575178403017 0.028574999999999999706 1 14 0.7448003145382178225 0.39165216296968274889 0.028574999999999999706 1 15 0.7048559114700959638 0.27543075527243382483");

    //ts.fromString("16 0.028575 1 0 0.64406291597264420101 1.7277722384035922865 0.028574999999999999706 1 1 0.54559659467979126912 0.70297292275039979792 0.028574999999999999706 1 2 0.47645297237633110754 0.470236504163515312 0.028574999999999999706 1 3 0.58269206446175303427 0.62763167874333680718 0.028574999999999999706 1 4 0.29422708869306196666 0.38221842046823018002 0.028574999999999999706 1 5 0.79098411064916840374 0.5791170519989594645 0.028574999999999999706 1 6 0.78840346724371590525 0.40737734901576322377 0.028574999999999999706 1 7 0.17467437750786263262 0.3111181200305964123 0.028574999999999999706 1 8 0.37199899511874784386 0.32953007291187219607 0.028574999999999999706 1 9 0.41355716716312423964 0.18650046168547351755 0.028574999999999999706 1 10 0.59538797421528366627 0.45979374468949513188 0.028574999999999999706 1 11 0.35409418077988646933 0.39453693718148868941 0.028574999999999999706 1 12 0.50398134482025969128 0.091368573658218640854 0.028574999999999999706 1 13 0.64234837420834745014 0.36754643809466358562 0.028574999999999999706 1 14 0.72603056597802351213 0.34710465158111769579 0.028574999999999999706 1 15 0.66898865339276369912 0.29005036493851243939 ");

    //ts.fromString("16 0.028575 1 0 0.73727069856753146837 1.6725624274018191517 0.028574999999999999706 1 1 0.68497569279095227301 0.58341889284069858856 0.028574999999999999706 1 2 0.39289258559348277933 0.33074303741881128094 0.028574999999999999706 1 3 0.33445123399281617704 0.43977414707707573882 0.028574999999999999706 1 4 0.43946715355428916938 0.43222968321841664929 0.028574999999999999706 1 5 0.74795674302126025079 0.41477637382461907123 0.028574999999999999706 1 6 0.60952984800542864274 0.41111767291977951988 0.028574999999999999706 1 7 0.37902485334936664163 0.39767547648851719888 0.028574999999999999706 1 8 0.49505597019927854996 0.44695789592834328507 0.028574999999999999706 1 9 0.49665153292490166193 0.31350063140152650831 0.028574999999999999706 1 10 0.55812384264767278541 0.12923725071253075458 0.028574999999999999706 1 11 0.11959690887785458246 0.19914007834707275757 0.028574999999999999706 1 12 0.3591223062960774981 0.15840918692849786931 0.028574999999999999706 1 13 0.61602960168335174629 0.46806477660225515036 0.028574999999999999706 1 14 0.65099658225470324879 0.18232981961394337245 0.028574999999999999706 1 15 0.69074021849157274389 0.085799138901048566108");

    //ts.fromString("16 0.028575 1 0 0.47754896338647262466 1.5327010250751615938 0.028574999999999999706 1 1 0.49023545940865820292 0.456176120793257589 0.028574999999999999706 1 2 0.36366093435322260907 0.52381734911634869611 0.028574999999999999706 1 3 0.69239686842708947356 0.45309238535245821566 0.028574999999999999706 1 4 0.33325791812688215776 0.39666323409204418127 0.028574999999999999706 1 5 0.64078361742171374704 0.27117231426935095584 0.028574999999999999706 1 6 0.77112116546391029637 0.35326092696568167639 0.028574999999999999706 1 7 0.49151363239183848419 0.28545493582668646226 0.028574999999999999706 1 8 0.34033400849553263612 0.25644485685156070964 0.028574999999999999706 1 9 0.64403103046620424621 0.38510597096073556633 0.028574999999999999706 1 10 0.59261162696567648567 0.41406185671869155662 0.028574999999999999706 1 11 0.4399556401873104261 0.3634794475349590881 0.028574999999999999706 1 12 0.23081752717133455222 0.19430661870844981309 0.028574999999999999706 1 13 0.59792469963788374976 0.10795690126859040781 0.028574999999999999706 1 14 0.70124955229238750132 0.39657676771855637776 0.028574999999999999706 1 15 0.87463099556027235071 0.38141830105100960058");

    //ts.fromString("16 0.028575 1 0 0.36599230818683831101 1.621480144385664568 0.028574999999999999706 1 1 0.5535410799184427022 0.7884059557385033612 0.028574999999999999706 1 2 0.48943818499774249808 0.50331109625704939514 0.028574999999999999706 1 3 0.6307555079877547044 0.6582416984032576357 0.028574999999999999706 1 4 0.55599670812713974932 0.56170436620102992542 0.028574999999999999706 1 5 0.47143832921793371593 0.41529130320927681863 0.028574999999999999706 1 6 0.72127535056373526245 0.55449866863692143237 0.028574999999999999706 1 7 0.60229010263725712981 0.37988092548559260209 0.028574999999999999706 1 8 0.54296558882115553146 0.304226498715007454 0.028574999999999999706 1 9 0.48160510901934816541 0.22336697686125550621 0.028574999999999999706 1 10 0.56667108102124863489 0.20512683136153220254 0.028574999999999999706 1 11 0.25806218173985251418 0.17892938889564119487 0.028574999999999999706 1 12 0.37492642398036740703 0.10161237902914954656 0.028574999999999999706 1 13 0.64274374239697285027 0.14723613200999760564 0.028574999999999999706 1 14 0.77851880120475647207 0.32922772504238806413 0.028574999999999999706 1 15 0.90002681314741406204 0.17413899704197483009");

    //ts.fromString("16 0.028575 1 0 0.36599230818683831101 1.621480144385664568 0.028574999999999999706 1 1 0.5535410799184427022 0.7884059557385033612 0.028574999999999999706 1 2 0.48943818499774249808 0.50331109625704939514 0.028574999999999999706 1 3 0.6307555079877547044 0.6582416984032576357 0.028574999999999999706 1 4 0.55599670812713974932 0.56170436620102992542 0.028574999999999999706 1 5 0.47143832921793371593 0.41529130320927681863 0.028574999999999999706 1 6 0.72127535056373526245 0.55449866863692143237 0.028574999999999999706 1 7 0.60229010263725712981 0.37988092548559260209 0.028574999999999999706 1 8 0.54296558882115553146 0.304226498715007454 0.028574999999999999706 1 9 0.48160510901934816541 0.22336697686125550621 0.028574999999999999706 1 10 0.56667108102124863489 0.20512683136153220254 0.028574999999999999706 1 11 0.25806218173985251418 0.17892938889564119487 0.028574999999999999706 1 12 0.37492642398036740703 0.10161237902914954656 0.028574999999999999706 1 13 0.64274374239697285027 0.14723613200999760564 0.028574999999999999706 1 14 0.77851880120475647207 0.32922772504238806413 0.028574999999999999706 1 15 0.90002681314741406204 0.17413899704197483009");

    //ts.fromString("16 0.028575 1 0 0.73886180524480260523 1.689938728050929706 0.028574999999999999706 1 1 0.6899839859204864867 0.55210595355856617594 0.028574999999999999706 1 2 0.54056860892922309336 0.4153083564464826738 0.028574999999999999706 1 3 0.62221070638522069096 0.62068998888005344483 0.028574999999999999706 1 4 0.4841816472557238793 0.32867069581445662596 0.028574999999999999706 1 5 0.43316084953394273782 0.45500222538141804618 0.028574999999999999706 1 6 0.58937579706155152071 0.52446986226645098128 0.028574999999999999706 1 7 0.30735196204698034261 0.25650917200768585946 0.028574999999999999706 1 8 0.37678423445265901703 0.26148106782759994449 0.028574999999999999706 1 9 0.65072860794990949884 0.17494046218490716393 0.028574999999999999706 1 10 0.84248841015524000841 0.31255724979970256339 0.028574999999999999706 1 11 0.4492503499209424267 0.54757214553059418982 0.028574999999999999706 1 12 0.44385149137846496847 0.18731003327754600085 0.028574999999999999706 1 13 0.33790889488430492227 0.15586721644424836319 0.028574999999999999706 1 14 0.5211641057416196432 0.50505054737842980384 0.028574999999999999706 1 15 0.7161328686306938085 0.20145561688845531334");

    //ts.fromString("16 0.028575 1 0 0.49172249518795058121 1.5665211194570352049 0.028574999999999999706 1 1 0.43609092112649311401 0.51490128696476278325 0.028574999999999999706 1 2 0.50605596078652437253 0.62109553529446182019 0.028574999999999999706 1 3 0.48087558719882128599 0.76964662455898869009 0.028574999999999999706 1 4 0.33458408738197892296 0.47246686717942842915 0.028574999999999999706 1 5 0.62036037567152491068 0.43747277117281607728 0.028574999999999999706 1 6 0.56285308750231677344 0.45194625833269480575 0.028574999999999999706 1 7 0.4985652031383188687 0.29504190087815829191 0.028574999999999999706 1 8 0.47052536272269923634 0.36214157901156784902 0.028574999999999999706 1 9 0.49177092003804540044 0.46326495142929186022 0.028574999999999999706 1 10 0.70724900429189596629 0.29746115237146814048 0.028574999999999999706 1 11 0.41489103903843937982 0.15015187488908801616 0.028574999999999999706 1 12 0.6291741637560859246 0.25442380956049026608 0.028574999999999999706 1 13 0.56800458695844402435 0.20285351943568244448 0.028574999999999999706 1 14 0.79693311101527175566 0.22464111837256789395 0.028574999999999999706 1 15 0.82602335091634804254 0.12847578678623400306");

    //ts.fromString("16 0.028575 1 0 0.36599230818683831101 1.621480144385664568 0.028574999999999999706 1 1 0.5535410799184427022 0.7884059557385033612 0.028574999999999999706 1 2 0.48943818499774249808 0.50331109625704939514 0.028574999999999999706 1 3 0.6307555079877547044 0.6582416984032576357 0.028574999999999999706 1 4 0.55599670812713974932 0.56170436620102992542 0.028574999999999999706 1 5 0.47143832921793371593 0.41529130320927681863 0.028574999999999999706 1 6 0.72127535056373526245 0.55449866863692143237 0.028574999999999999706 1 7 0.60229010263725712981 0.37988092548559260209 0.028574999999999999706 1 8 0.54296558882115553146 0.304226498715007454 0.028574999999999999706 1 9 0.48160510901934816541 0.22336697686125550621 0.028574999999999999706 1 10 0.56667108102124863489 0.20512683136153220254 0.028574999999999999706 1 11 0.25806218173985251418 0.17892938889564119487 0.028574999999999999706 1 12 0.37492642398036740703 0.10161237902914954656 0.028574999999999999706 1 13 0.64274374239697285027 0.14723613200999760564 0.028574999999999999706 1 14 0.77851880120475647207 0.32922772504238806413 0.028574999999999999706 1 15 0.90002681314741406204 0.17413899704197483009"); 

    //ts.fromString("16 0.028575 1 0 0.60099903563179357668 1.6052197366602589668 0.028574999999999999706 1 1 0.54395965425949144301 0.59114385647648126643 0.028574999999999999706 1 2 0.40163049078661816615 0.61821528143984805226 0.028574999999999999706 1 3 0.61284048199552754177 0.46597269212189545984 0.028574999999999999706 1 4 0.39591436443006311485 0.40483649596417220495 0.028574999999999999706 1 5 0.47195722804196149625 0.34220761125971738137 0.028574999999999999706 1 6 0.78960830059573539064 0.46090371313268935216 0.028574999999999999706 1 7 0.30476673030854289914 0.37331541464864820279 0.028574999999999999706 1 8 0.59399134364894712323 0.36369158561944114894 0.028574999999999999706 1 9 0.54026130594422250297 0.40430633467509852208 0.028574999999999999706 1 10 0.76340358951142039956 0.36610586752945561972 0.028574999999999999706 1 11 0.35026345873097419759 0.29959253678500574747 0.028574999999999999706 1 12 0.47758197915555644641 0.42350833437238205592 0.028574999999999999706 1 13 0.65380290142999886172 0.13891442781190763145 0.028574999999999999706 1 14 0.71275194898540950028 0.3083551919261416363 0.028574999999999999706 1 15 0.69098056400282170664 0.22890993460364095213"); 

    //ts.fromString("16 0.028575 1 0 0.45625642400992555414 1.7456630140759121783 0.028574999999999999706 1 1 0.49842646103815851921 0.59962010592250414298 0.028574999999999999706 1 2 0.35725848934570525461 0.45466688347352995914 0.028574999999999999706 1 3 0.68318229226172400015 0.43193182529731666275 0.028574999999999999706 1 4 0.45444510014735406411 0.41714923549542287651 0.028574999999999999706 1 5 0.56854711853717898595 0.57725355423796820276 0.028574999999999999706 1 6 0.5203120995044332453 0.4257720606446704914 0.028574999999999999706 1 7 0.48984383826683008945 0.23935063208251985967 0.028574999999999999706 1 8 0.59334545016570827691 0.4659455965825816115 0.028574999999999999706 1 9 0.60794742497407083803 0.27979313656583204573 0.028574999999999999706 1 10 0.77210914703640765033 0.60526130976043568399 0.028574999999999999706 1 11 0.42690320087955191397 0.24114013157354877159 0.028574999999999999706 1 12 0.27623501589937560219 0.13692707019774241761 0.028574999999999999706 1 13 0.58000392041637982565 0.38697640224404022957 0.028574999999999999706 1 14 0.71593041399601053953 0.089190085190408871507 0.028574999999999999706 1 15 0.82869056052852407834 0.24873547155052042057");
    //cerr << "rack:" << ts.isValidBallPlacement() << endl;
}

void printEvents(list<Event*> events) {
    for (list<Event*>::iterator e = events.begin(); e != events.end(); ++e) {
        cerr << (**e) << endl;
    }
}

void printOverlap(TableState& ts) {
    for (vector<Ball>::iterator b1 = ts.getBegin(); b1 != ts.getEnd()-1; ++b1) {
        for (vector<Ball>::iterator b2 = b1+1; b2 != ts.getEnd(); ++b2) {
            double ovl=(b1->getRadius()+b2->getRadius())-
                sqrt(pow(b1->getPos().x - b2->getPos().x, 2.0) + pow(b1->getPos().y - b2->getPos().y,2.0));
            if (ovl>0) {
                cerr << b1->getIDString() << " overlaps " << b2->getIDString() << endl;
                cerr << "overlap degree: " << ovl << endl;
            }
        }
    }
}

TableState* getTestState() {
    TableState* ts = new TableState();
    rack(*ts);
    return ts; 
}

ShotParams* getTestShotParams() {
    ShotParams* sp = new ShotParams(0,0,5,270,5);
    return sp;
}

gsl_rng *Utils::_rng= NULL;

gsl_rng * Utils::rng()
{
    if (!_rng) {
        gsl_rng_env_setup();
        const gsl_rng_type * rngType = gsl_rng_default;
        _rng = gsl_rng_alloc(rngType);
        if (!gsl_rng_default_seed) {
            gsl_rng_set(_rng, time(0));
        }
    }
    return _rng;
}

void Ball::moveAway(Ball& other, double epsilon)
{
    r=r+((r-other.r).to_v().norm()*epsilon).to_p();
}


void TableState::fixOverlap(const bool VERBOSE)
{
    int count=100;
    bool overlap=false;
    do {
        overlap=false;
        for(vector<Ball>::iterator i = getBegin(); i != getEnd(); ++i)
            if (i->isInPlay())
            {
                for(vector<Ball>::iterator j = getBegin(); j != i; ++j)
                    if (j->isInPlay()) {
                        //int moveCount=100;
                        if (i->overlaps(*j,0)/* && --moveCount*/) {
                            overlap=true;
                            if (VERBOSE) {
                                printOverlap(*this);
                                cerr << "Moving Back balls (" << count << ") " << *i << *j << endl;
                            }
                            i->moveAway(*j,Utils::EPSILON);
                            j->moveAway(*i,Utils::EPSILON);
                            if (VERBOSE) {
                                cerr << "After moveback " << *this;
                            }
                        }
                    }
            }
        //if (count<99) cerr << '.';
    } while (overlap && --count);
#ifdef POOLFIZ_ERRORS        
    if (!count) {
        for(vector<Ball>::iterator i = getBegin(); i != getEnd(); ++i)
            if (i->isInPlay())
            {
                for(vector<Ball>::iterator j = getBegin(); j != i; ++j)
                    if (j->isInPlay())
                    {
                        if (i->overlaps(*j,0)) {
                            printOverlap(*this);
                            cerr << *i << " still overlaps " << *j << endl;
                            toStream(cerr);
                            POOLFIZ_ERROR();
                        }
                    }
            }
    }
#endif        
}


}//namespace pool
