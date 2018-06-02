#ifndef INFOBOX_HH
#define INFOBOX_HH

// C++ includes
#include <iostream>
#include <string>

// ROOT includes
#include "TPaveText.h"
#include "TText.h"
#include "TString.h"

/** @brief Stores the MICE run information to be printed as TPaveText.
 *
 *	   Defines the style of the MICE information text box. Sets the default position of the
 *	   information box, which may be overwritten if necessary.
 */
class InfoBox {
 public:

  /** @brief MAUS info box constructor, sets the necessary information
   *
   *  @param	type		Type of data presented: [simulated], Preliminary, etc.
   *  @param	maus		Version of MAUS used
   *  @param	run		Run simulated or reconstructed
   *  @param	cycle		ISIS user cycle tag
   */
  InfoBox(const std::string type="",
	  const std::string maus="",
	  const std::string run="",
	  const std::string cycle="");

  /** @brief Copy constructor */
  InfoBox(const InfoBox& ib);

  /** @brief Equality operator */
  InfoBox& operator=(const InfoBox& ib);

  /** @brief Destructor */
  ~InfoBox();

  /** @brief Updates the position of the info box
   *
   *  @param	pos		Corner in which to draw the info box (tl, tr, bl, br)
   *  @param	alpha		Relative margin from the edges of the TCanvas [%]
   *  @param	beta		Relative size [%]
   */
  void SetPosition(const std::string& pos,
		   const double alpha=.02,
		   const double beta=.16);

  /** @brief Updates the position of the info box
   *
   *  @param	xmin		Minimum x position in NDC
   *  @param	ymin		Minimum y position in NDC
   *  @param	xmax		Maximum x position in NDC
   *  @param	ymax		Maximum y position in NDC
   */
  void SetPosition(const double& xmin,
		   const double& ymin,
		   const double& xmax,
		   const double& ymax);

  /** @brief Sets the opacity of the TPaveText (0 is transparent, 1 is fully opaque)
   *
   *  @param	opacity		Value of the opacity (from 0 to 1)
   */
  void SetOpacity(const double& opacity);

  /** @brief Draws the info box on the current canvas */
  void Draw();

 private:

  TPaveText*		_text;		///< Info box
  std::string		_type;		///< Type of data represented
  std::string		_maus;		///< Version of MAUS
  std::string		_run;		///< Run tag
  std::string		_cycle;		///< ISIS user cycle
};

#endif
