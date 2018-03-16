/* AIBase.h
 * -------------------
 * This file is the header for the AI's
 * Every AI should inherit from this abstract class.
 */


#ifndef _AIBASE_H_
#define _AIBASE_H_


#include <stdio.h>
#include <string>
#include <libconfig.hh>
#include "Rules.h"
#include "Noise.h"
#include "Stopwatch.h"

using namespace std;

#define MAX_RETRIES 10

namespace Pool {

    /** Base class for pool AI implementations.
     *  When writing your AI, inherit from this class.
     */
    class AIBase
    {
        public:
            /** Constructor. Generates a new AI instance.
             *  \param config Config object with configuration information for the AI.
             *  \param stopwatch A Stopwatch object allowing access to the time model used.
             */
            AIBase(libconfig::Config &config, Stopwatch* stopwatch): _config(config), _stopwatch(stopwatch) {};

            /** Compute a shot for a given game state and noise information.
             *  Default implementation keeps shot information, makes sure the AI can handle the game type, 
             *  and then calls decide() or shoot() based on the turn type.
             * 
             *  Generally, you need not re-implement this method in your AI.
             */
            virtual GameShot& computeShot(const GameState &gs, const Pool::Noise* const noise=NULL);

            /** Recompute a shot. Called after a shot sent by the AI was deemed physcially impossible or underspecified.
             *  Default implementation retries previous shot again up to MAX_RETRIES times, then concedes the game.
             * 
             *  Generally, you need not re-implement this method in your AI.
             */
            virtual GameShot& reComputeShot();

            /** Return the name of your AI.
             *  Default implementation reads it from the config file under "name".
             * 
             *  Generally, you need not re-implement this method in your AI.
             */
            virtual string getName() const {return _config.lookup("name");};

            /** Return the password for your AI.
             *  Default implementation reads it from the config file under "password".
             * 
             *  Generally, you need not re-implement this method in your AI.
             */
            virtual string getPassword() const  {return _config.lookup("password");};

            /** Return the owner ID for the AI.
             *  Owner ID is used in a multi-user server to allow non-admin users to view information related to the agent.
             *  Default implementation does not specify an owner ID.
             * 
             *  Generally, you need not re-implement this method in your AI.
             */
            virtual string getOwner() const {return "";};

            /** Check if the AI can handle a specific game type.
             *  Default implementation is to always return false.
             *  
             *  AI authors should re-impelment this method to return true for game types that are supported by the AI.
             */
            virtual bool forGame(GameType gt) {return false;}

        protected:
            /** Make a decision when a decision node is reached. State and noise information is stored in class memebers.
             * Default implementation is to rerack whenever an eight-ball decision is needeed.
             * 
             * AI authors should re-implement this method if any other decision is desired.
             */
            virtual Decision decide();
            /** Make a shot when one is needed.
             * The game state and noise information are stored in class members. Resulting shot should be written to the "shot" member.
             * Default implementation calls breakShot() or otherShot() depending on the turn type.
             * 
             * Generally, you need not re-implement this method in your AI.
             */
            virtual void shoot();

            /** Make a break shot.
             * This is a pure virtual function that must be re-implemented by any AI. 
             * The game state and noise information are stored in class members. Resulting shot should be written to the "shot" member.
             */
            virtual void breakShot() = 0;
            /** Make a non-break shot. This includes ball-in-hand shots.
             * This is a pure virtual function that must be re-implemented by any AI. 
             * The game state and noise information are stored in class members. Resulting shot should be written to the "shot" member.
             */
            virtual void otherShot() = 0;

            /** Configuration information for the agent */
            libconfig::Config& _config;

            /** Stopwatch object defining time model used by the agent */
            Stopwatch* const _stopwatch;

            /** State of the game when AI is needed */
            const GameState *gameState;

            /** Noise information */
            const Noise *noise;

            /** Shot computed by AI */
            GameShot shot;

            /** Number of attempts remaining for an illegal/impossible shot */
            int retries;
    };

};

#endif
