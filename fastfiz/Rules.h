#ifndef POOL_RULES_H
#define POOL_RULES_H

#include "FastFiz.h"
#ifdef SWIG
%include "std_string.i"
%{
#include "Rules.h"
    using namespace Pool;
    %}
#endif /* SWIG */
#include <string>
    using namespace std;

    namespace Pool {
        /**
         * A definition of all the possible decisions a client can make.
         * If a state does not require a decision, but rather only a shot,
         * set the decision of a GameShot to DEC_NO_DECISION. Otherwise, use
         * the appropriate decision.
         */
        enum Decision {
            DEC_NO_DECISION, /**< Default decision when none is needed */
            DEC_KEEP_SHOOTING, /**< Keep your turn after opponent's bad break or push out */
            DEC_RERACK, /**< Re-rack the balls and get a turn after opponent's bad break */
            DEC_EIGHTBALL_RERACK_OPP_SHOOT, /**< If TT_EIGHTBALL_FOUL_ON_BREAK, the incoming player can decide to re-rack and have the opponent (original breaker) re-break. If they wimply want to re-rack and themselves shoot, choose DEC_RERACK*/
            DEC_CONCEDE, /**< concede the game, always allowed */
            DEC_NINEBALL_PUSH_OUT /**< In TT_NINEBALL_FIRST_SHOT, means declaring "Push Out". In TT_NINEBALL_PUSH_OUT means player does not accept position. */
        };

        /**
         * An enumeration of all the possible turn types an agent might have to 
         * deal with. The turn types are as general as possible, with a few additions
         * as required by each specific game.
         */
        // Turn Types -- what kind of shot does the player have
        enum TurnType {
            TT_NORMAL, /**< Standard shot, cue ball must not be placed */
            TT_BALL_IN_HAND, /**< Cue ball must be placed anywhere on the table */
            TT_BEHIND_LINE, /**< Cue ball must be placed behind the headstring */
            TT_RESERVED, /**< Not used */
            TT_BREAK, /**< Break shot. Cue ball must be placed behind the headstring. */
            TT_WIN, /**< Current player wins. This is a terminal state. */
            TT_EIGHTBALL_FOUL_ON_BREAK, /**< If there is a foul on the break. Needs a decision */
            TT_EIGHTBALL_8BALL_POCKETED_ON_BREAK, /**< If the eight-ball is pocketed on the break. Needs a decision */
            TT_NINEBALL_FIRST_SHOT, /**< First shot after a legal break in nine-ball */
            TT_NINEBALL_PUSH_OUT /**< Incoming player's decision node after Push Out (neesds to decide DEC_KEEP_SHOOTING or DEC_NINEBALL_PUSH_OUT) */
        };

        inline istream & operator>> (istream & is, TurnType& tt)
        {
            int code;
            is >> code;
            tt = TurnType(code);
            return is;
        }

        // Game Types -- what billards game are we playing
        /**
         * An enumeration of all the game types. Only Eight-Ball is supported.
         */
        enum GameType {GT_NONE,GT_EIGHTBALL,GT_NINEBALL,GT_SNOOKER,GT_ONEPOCKET};

        // Shot Results -- any pool shot that is actually executed is "OK", 
        //                 even if leads to loss of game.
        /**
         * The set of all possible results of an execution of a shot. 
         */
        enum ShotResult {
            SR_OK, /**< The shot was executed, and the shooter won the game or may shoot again. */
            SR_OK_LOST_TURN, /**< The shot was executed, but the shooter has lost their turn, and possibly the game. */
            SR_BAD_PARAMS, /**< The shot information given is not applicable to the game state. */
            SR_SHOT_IMPOSSIBLE, /**< The shot is not physically possible. For example, if a ball is in the way of the cue stick. */
            SR_TIMEOUT /**< The agent took too much time computing the shot, and timed out. Shot was not executed. */
        };

        /**
         * The object that holds all parameters necessary for a complete shot.
         */  

        struct GameShot {
            ShotParams params; /**< Struct of the shot parameters */
            double cue_x,cue_y; /**< Cue ball positioning for ball-in-hand */
            Ball::Type ball; /**< Called ball */
            Table::Pocket pocket; /**< Called pocket */
            Decision decision; /**< The decision (usually not required) */
            /** Amount of time spent by the AI computing the shot. This value is ignored during competition
             * however if the "fake time" option is set on the server, than the server uses
             * this value. (This is used for testing environments, especially when one
             * does not have full access to an unloaded computer). */
            double timeSpent;
        };

        /**
         * This defines most of the information and functions necessary about the state
         * of a game. All rules objects subclass off this.
         */
        // Abstract base class, signifies all information about the state of a game.
        class GameState {
            public:
                ////////////////////////////////////////////////////////
                // FACTORY FUNCTION:
                // This is the appropriate way to construct the GameStates
                // From a source stream
                ////////////////////////////////////////////////////////
                /**
                 * Generate a new GameState object by reading information from an input stream.
                 * Game state should have been written with toStream or toString.
                 */
                static GameState* Factory(istream & sourceStream);
                /**
                 * Generate a new GameState object by reading information from an input string.
                 * Game state should have been written with toStream or toString.
                 */
                static GameState* Factory(string sourceString)
                {
                    istringstream gameSourceStream(sourceString);
                    return Factory(gameSourceStream);
                }
                /**
                 * Returns a newly allocated racked (initial) state of a game by game type.
                 */
                static GameState* RackedState(GameType gameType);


#ifndef SWIG


                //////////////////////////////////////////////////////////
                // CONSTRUCTORS
                // default copy constructor ok
                //////////////////////////////////////////////////////////
                /** Create new game state. 
                 * This function should not be called directly. Use Factory instead.
                 * Zero timeouts (default) mean no timeout. 
                 */
                explicit GameState (double timeleft =0.0, double timeleft_opp=0.0) : 
                    _tableState(), _turnType(TT_BREAK), 
                    _timeLeft (timeleft), _timeLeftOpp (timeleft_opp), 
                    _curPlayerStarted(true) 
            {}; 

                /** Parse common part of GameState from string */
                GameState (string gameString)
                {
                    istringstream gameStream(gameString);
                    importFromStream(gameStream);
                };

                /** Default destructor */
               virtual ~GameState() { };
                
                /**
                 * Determines if the GameState is in a terminal state (i.e. the game is over.)
                 */
                virtual bool isTerminal() const {return _turnType == TT_RESERVED || _turnType == TT_WIN;};
                /**
                 * Determines if the agent must provide a valid position for the cue ball.
                 */
                virtual bool positionRequired() const {return _turnType == TT_BREAK || _turnType == TT_BALL_IN_HAND || _turnType == TT_BEHIND_LINE;};
                /**
                 * Determines if the agent is allowed to make a decision.
                 */
                virtual bool decisionAllowed() const {return !shotRequired() && !isTerminal();};
                /**
                 * Determines if the agent is required to provide a valid shot.
                 */
                virtual bool shotRequired() const {return !isTerminal();};

                ///////////////////////////////////////////////////
                // Input and output funzies.
                ///////////////////////////////////////////////////
                // Define friend operator functions
                /**
                 * Convienence parsing operators.
                 */
                friend ostream& operator<<(ostream& os, const GameState& obj);
                friend istream& operator>>(istream& is, GameState& obj);
#endif /* ! SWIG */

                // Serilaize GameState as a string, only relevant parts.
                // Note: base classes can override, but default should be good
                /**
                 * Serializes the GameState as a string.
                 */
                virtual string toString();

#ifndef SWIG
                // Puts out the default values.
                // Subclasses should call this, and then output their internals
                /**
                 * Exports the GameState to a stream.
                 */
                virtual void toStream(ostream & out) const ;

                /** Allocate a new copy of this GameState */
                virtual GameState* clone () = 0;

#endif /* ! SWIG */
                // Return type of game
                // Should be overridden in lower classes
                /**
                 * Returns the GameType of the GameState.
                 */
                virtual GameType gameType() const  =0;


                //////////////////////////////////////////////
                // AI FUNCTIONS
                /////////////////////////////////////////////

                ////////////////////////////////////////////
                // Get various properties of the game state
                ///////////////////////////////////////////
                /** Returns true if solids/stripes have not been determined yet */
                virtual bool isOpenTable() const {return true;};
                virtual TurnType getTurnType() const {return _turnType;};
                /** Returns true if current player should shoot solids. 
                 * Only meaningful if isOpenTable is false. 
                 */
                virtual bool playingSolids() const {return false;};
                /** Returns true if the current player had the inital break */
                virtual bool curPlayerStarted() const {return _curPlayerStarted;};
                /** Returns the amount of time left (in seconds) for the current player */
                virtual double timeLeft() const {return _timeLeft;};
                /** Returns the amount of time left (in seconds) for the other player */
                virtual double timeLeftOpp() const {return _timeLeftOpp;};
                /** Allows caller read-only access to the internal table state. */
                virtual const TableState& tableState() const {return _tableState;};


                // Get Result of shot. In base class, this checks for timeout, calls isAllowed, then tries
                // to execute the shot with the physics, then finally calls processEventList.
                /**
                 * Calls the rules to check the shot values, executes the shot using the 
                 * physics library, and then evalues the resulting GameState.
                 * Optionally, returns the Shot object used to execute the physical shot.
                 * \param shot Shot to execute
                 * \param shotObj Pointer to an unallocated Shot* which will include a newly allocated Shot object.
                 *                Caller is responsible to delete *shotObj after processing.
                 */
                virtual ShotResult executeShot(const GameShot& shot, Shot **shotObj=NULL);      

#ifndef SWIG


            protected:

                /** Returns true if the current turnType is defined for the current game. */
                virtual bool isLegalTurnType() const {
                    return _turnType==TT_NORMAL || _turnType==TT_BALL_IN_HAND || 
                        _turnType == TT_BEHIND_LINE || _turnType==TT_BREAK  || 
                        _turnType==TT_WIN;
                }

                /** Read common information from an input stream */
                virtual void importFromStream(istream & sourceStream);

                //////////////////////////////////////////////////
                // Rules implementation functions
                //////////////////////////////////////////////////
                /**
                 * Possible results of pre-evaluation of shot parameters for current game state.
                 */
                enum PreProcessCode { 
                    PPC_BADPARAMS=0, /**< Shot parameters not satisfactory for game state. */
                    PPC_G_NOEXECUTE, /**< Decision information in shot is sufficient. No need to execute shot. */
                    PPC_G_NORMAL,    /**< Shot should be processed normally. Cue ball is not to be moved. */
                    PPC_G_PLACECUE   /**< Shot should be processed normally. Cue ball is to be placed at requested position. */
                };
                /** Pre-evaluate shot parameters for the current game state */
                virtual PreProcessCode preProcess(const GameShot& shot) =0;

                /** Called after shot execution to process an event list and adjust game state. */
                virtual void processShot(const vector<Event*>& eventList, const GameShot & gameShot) =0;

                /** Update TableState to the racked state, does NOT update other params. */
                virtual void rack()=0;

                /** Utility function to switch the active player */
                virtual void switchSides();

                //////////////////////////////////////////////
                // INTERNAL VARIABLES
                //////////////////////////////////////////////
                TableState _tableState;
                TurnType _turnType;
                double _timeLeft,_timeLeftOpp; // in seconds
                bool _curPlayerStarted;
                bool _switchedSides;

#endif /* ! SWIG */
        }; // Class GameState

#ifndef SWIG

        /** Implementation of the rules of Eight-Ball */
        class EightBallState : public GameState
        {
            public:

                /** Create new (racked) Eight Ball state.*/
                explicit EightBallState (double timeleft =0.0, double timeleft_opp=0.0) : 
                    GameState(timeleft,timeleft_opp), _openTable(true), _solids(true)
            {
                rack();
            };


                // The usual IO friend functions
                /**
                 * Convienence functions to facilitate working with the strings used to 
                 * send over the wire.
                 */
                friend ostream& operator<<(ostream& os, const EightBallState& obj);
                /**
                 * See above
                 */
                friend istream& operator>>(istream& is, EightBallState& obj);

                /**
                 * Construct an EightBallState by reading from string. Used by Factory.
                 */
                EightBallState (string gameString)
                {
                    istringstream infoStream(gameString);
                    importFromStream(infoStream);
                };

                // Return type of game
                /**
                 * Returns the type of game (GT_EIGHTBALL)
                 */
                GameType gameType() const  {return GT_EIGHTBALL;};

                // Parse game specific parts of the stream
                /**
                 * Construct an EightBallState by reading from stream. Used by Factory.
                 */
                EightBallState (istream &is) { importFromStream(is);};

#ifndef SWIG
                // Puts out the default values.
                // Subclasses should call this, and then output their internals
                virtual void toStream(ostream & out) const ;

                virtual GameState* clone ();

#endif /* ! SWIG */

                virtual bool playingSolids() const {return _solids;};
                virtual bool isOpenTable() const {return _openTable;};

                virtual bool shotRequired() const {
                    return GameState::shotRequired() && _turnType!=TT_EIGHTBALL_8BALL_POCKETED_ON_BREAK &&
                        _turnType!=TT_EIGHTBALL_FOUL_ON_BREAK;
                };


            protected:

                // Import from a stream
                virtual void importFromStream(istream & sourceStream);

                virtual bool isLegalTurnType() const {
                    return GameState::isLegalTurnType() || 
                        _turnType==TT_EIGHTBALL_8BALL_POCKETED_ON_BREAK ||
                        _turnType==TT_EIGHTBALL_FOUL_ON_BREAK;
                }


                ////////////////////////////////////////////////////////
                // The actual rules implementation functions
                ////////////////////////////////////////////////////////

                // Check if given shot supplies all required parameters for current turnType.
                // Returns an emum: bad params, don't execute but good, and execute, execute + cue position
                virtual PreProcessCode preProcess(const GameShot& shot);

                // Update turnType and curPlayerStarted based on the events that occured in the shot
                virtual void processShot(const vector<Event*>& eventList, const GameShot & gameShot);

                // Update TableState to the racked state, does NOT update other params.
                virtual void rack();

                virtual void switchSides();

                /////////////////////////////////////////////////////////////////////
                // Internal Variables
                /////////////////////////////////////////////////////////////////////
                bool _openTable, _solids;


        }; // Class EightBallState

        /** Implementation of the rules of Nine-Ball */
        class NineBallState : public GameState
        {
            public:

                // Create new game state. Zero timeouts (default) mean no timeout.
                explicit NineBallState (double timeleft =0.0, double timeleft_opp=0.0) : 
                    GameState(timeleft,timeleft_opp), _fouls(0), _opp_fouls(0)
            {
                rack();
            };


                // The usual IO friend functions
                friend ostream& operator<<(ostream& os, const NineBallState& obj);
                friend istream& operator>>(istream& is, NineBallState& obj);

                // Parse common part of GameState from string
                NineBallState (string gameString)
                {
                    istringstream infoStream(gameString);
                    importFromStream(infoStream);
                };

                // Return type of game
                /**
                 * Returns the type of game (GT_NINEBALL)
                 */
                GameType gameType() const  {return GT_NINEBALL;};

                // Parse game specific parts of the stream
                NineBallState (istream &is) { importFromStream(is);};

#ifndef SWIG
                // Puts out the default values.
                // Subclasses should call this, and then output their internals
                virtual void toStream(ostream & out) const ;

                virtual GameState* clone ();

#endif /* ! SWIG */

                virtual bool shotRequired() const {
                    return GameState::shotRequired() && _turnType != TT_NINEBALL_PUSH_OUT;
                };


            protected:

                // Import from a stream
                virtual void importFromStream(istream & sourceStream);

                virtual bool isLegalTurnType() const {
                    return _turnType != TT_BEHIND_LINE && (GameState::isLegalTurnType() || _turnType == TT_NINEBALL_FIRST_SHOT || _turnType == TT_NINEBALL_PUSH_OUT); 
                }


                ////////////////////////////////////////////////////////
                // The actual rules implementation functions
                ////////////////////////////////////////////////////////

                // Check if given shot supplies all required parameters for current turnType.
                // Returns an emum: bad params, don't execute but good, and execute, execute + cue position
                virtual PreProcessCode preProcess(const GameShot& shot);

                // Update turnType and curPlayerStarted based on the events that occured in the shot
                virtual void processShot(const vector<Event*>& eventList, const GameShot & gameShot);

                // Update TableState to the racked state, does NOT update other params.
                virtual void rack();

                virtual void switchSides();

                /////////////////////////////////////////////////////////////////////
                // Internal Variables
                /////////////////////////////////////////////////////////////////////
                int _fouls,_opp_fouls;


        }; // Class NineBallState

#endif /* ! SWIG */
        string getRulesVersion();

    } // Namespace pool

#endif
