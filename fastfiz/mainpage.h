/** 
  @mainpage Computational Billiards Library
  @author Alon Altman
  \section Introduction
  The FastFiz library consists of two major components: The physics engine and the rules engine, and several auxiliary components. The library also contains perl and python bindings for these components. The library also contains a base class for defining pool AIs (Pool::AIBase).

  \section FastFiz
  FastFiz.h implements the motion and collision physics of balls on a pool table independently of rules of any specific game. The Pool::TableState class represents a specific static position of any number of balls on a pool table. ShotParams is a simple struct containing the physical information regarding a cue strike. To simulate a shot, use Pool::TableState::executeShot on an existing Pool::TableState object. This function will modify the Pool::TableState to the final position of the shot and return newly allocated Shot object with the details of the shot. The Pool::Shot::getEventList function returns a list of Pool::Event objects which represent the various events that have occurred during the shot. These events are represented by the Pool::Event class hierarchy.

  \section Rules
  The rules library (Rules.h) implements the rules of billiards games by implementing an abstract Pool::GameState class capturing the positioning of the balls on the table and other information needed in the context of a billiards game. The main function used is Pool::GameState::executeShot (which expects a Pool::GameShot object) and modifies the Pool::GameState to the state after the shot has been executed. The returned Pool::ShotResult (in conjunction with examining the update Pool::GameShot object) provides information regarding the outcome of the shot.

  \section addcomp Additional components
  The Noise.h file defines a class hierarchy (Pool::Noise) for defining types of noise. The Stopwatch.h file defines a Stopwatch class hierarchy used for measuring time (real, CPU, or virtual) in the client and AI.

  The AiBase.h file defines Pool::AIBase, an abstract base class for Pool AIs. For a valid AI, you must minimally define forGame, breakShot, otherShot, and a constructor.

*/
