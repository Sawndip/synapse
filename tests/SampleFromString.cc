// C++ includes
#include "iostream"
#include "vector"
#include "algorithm"

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TRandom3.h"
#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoBBox.h"
#include "TGeoTube.h"
#include "TGeoArb8.h"
#include "TGeoMatrix.h"
#include "TView.h"

TGeoVolume* Letter(std::string& letter) {

  // Parent volume that contains all the elements  
  TGeoVolume* solid(NULL);
  TGeoBBox *box = new TGeoBBox(.5, .5, .1);
  solid = new TGeoVolume(letter.c_str(), box);
  solid->SetVisibility(kFALSE);

  if ( letter == "A" ) {
    // Add the two diagonals of the letter first
    TGeoArb8 *dbox = new TGeoArb8(.1);
    dbox->SetVertex(0, -.3, -.5);
    dbox->SetVertex(1, -.5, -.5);
    dbox->SetVertex(2, -.25, .5);
    dbox->SetVertex(3, -0.05, .5);
    dbox->SetVertex(4, -.3, -.5);
    dbox->SetVertex(5, -.5, -.5);
    dbox->SetVertex(6, -.25, .5);
    dbox->SetVertex(7, -0.05, .5);
    TGeoVolume *dbar = new TGeoVolume("diag", dbox);
    solid->AddNode(dbar, 1);
    TGeoRotation *rot = new TGeoRotation();
    rot->RotateY(180.);
    solid->AddNode(dbar, 2, new TGeoCombiTrans(-.3, 0., 0., rot));

    // Add two horizontal bar
    TGeoArb8 *hbox = new TGeoArb8(.1);
    hbox->SetVertex(0, -.05, -.3);
    hbox->SetVertex(1, -.25, -.3);
    hbox->SetVertex(2, -.2, -.1);
    hbox->SetVertex(3, -.1, -.1);
    hbox->SetVertex(4, -.05, -.3);
    hbox->SetVertex(5, -.25, -.3);
    hbox->SetVertex(6, -.2, -.1);
    hbox->SetVertex(7, -.1, -.1);
    TGeoVolume *hbar = new TGeoVolume("hor", hbox);
    solid->AddNode(hbar, 1);

  } else if ( letter == "B" ) {
    // Add the vertical bar of the letter
    TGeoBBox *vbox = new TGeoBBox(.1, .5, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.4, 0, 0.));

    // Add the three stubs to the right of the vertical bar
    TGeoBBox *sbox = new TGeoBBox(.1, .1, .1);
    TGeoVolume *sbar = new TGeoVolume("sbar", sbox);
    solid->AddNode(sbar, 1, new TGeoTranslation(-.2, 0., 0.));
    solid->AddNode(sbar, 2, new TGeoTranslation(-.2, .4, 0.));
    solid->AddNode(sbar, 3, new TGeoTranslation(-.2, -.4, 0.));

    // Add the two round parts
    TGeoTubeSeg *tub = new TGeoTubeSeg(.1, .3, .1, -90., 90.);
    TGeoVolume *round = new TGeoVolume("round", tub);
    solid->AddNode(round, 1, new TGeoTranslation(-.1, .2, 0.));
    solid->AddNode(round, 2, new TGeoTranslation(-.1, -.2, 0.));

  } else if ( letter == "C" ) {
    // Add two vertical bar to the side of the letter
    TGeoBBox *vbox = new TGeoBBox(.1, .15, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.4, 0., 0.));

    // Add the two round parts to the top and bottom
    TGeoTubeSeg *tub = new TGeoTubeSeg(.15, .35, .1, 0., 180.);
    TGeoVolume *round = new TGeoVolume("round", tub);
    solid->AddNode(round, 1, new TGeoTranslation(-.15, .15, 0.));
    TGeoRotation *rot = new TGeoRotation("rot", 90., 0., 90., 270., 0., 0.);
    solid->AddNode(round, 2, new TGeoCombiTrans(-.15, -.15, 0., rot));

  } else if ( letter == "D" ) {
    // Add a short vertical bar to the sides of the letter
    TGeoBBox *vrbox = new TGeoBBox(.1, .5, .1);
    TGeoVolume *vrbar = new TGeoVolume("vrbar", vrbox);
    solid->AddNode(vrbar, 1, new TGeoTranslation(-.4, 0., 0.));

    // Add the two stubs to the right of the vertical bar
    TGeoBBox *sbox = new TGeoBBox(.1, .1, .1);
    TGeoVolume *sbar = new TGeoVolume("sbar", sbox);
    solid->AddNode(sbar, 1, new TGeoTranslation(-.2, .4, 0.));
    solid->AddNode(sbar, 2, new TGeoTranslation(-.2, -.4, 0.));

    // Add a short vertical bar to the right of the letter
    TGeoBBox *vlbox = new TGeoBBox(.1, .15, .1);
    TGeoVolume *vlbar = new TGeoVolume("vlbar", vlbox);
    solid->AddNode(vlbar, 2, new TGeoTranslation(.1, 0., 0.));

    // Add the two round parts to the top and bottom
    TGeoTubeSeg *tub = new TGeoTubeSeg(.15, .35, .1, 0., 90.);
    TGeoVolume *round = new TGeoVolume("round", tub);
    solid->AddNode(round, 1, new TGeoTranslation(-.15, .15, 0.));
    TGeoRotation *rot = new TGeoRotation("rot", 90., 0., 90., 270., 0., 0.);
    solid->AddNode(round, 2, new TGeoCombiTrans(-.15, -.15, 0., rot));

  } else if ( letter == "E" ) {
    // Add the vertical bar of the letter
    TGeoBBox *vbox = new TGeoBBox(.1, .5, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.4, 0., 0.));

    // Add the horizontal bars of the letter
    TGeoBBox *hbox = new TGeoBBox(.25, .1, .1);
    TGeoVolume *hbar = new TGeoVolume("hbar", hbox);
    solid->AddNode(hbar, 1, new TGeoTranslation(-.05, .4, 0.));
    solid->AddNode(hbar, 2, new TGeoTranslation(-.05, 0., 0.));
    solid->AddNode(hbar, 3, new TGeoTranslation(-.05, -.4, 0.));

  } else if ( letter == "F" ) {
    // Add the vertical bar of the letter
    TGeoBBox *vbox = new TGeoBBox(.1, .5, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.4, 0., 0.));

    // Add the horizontal bars of the letter
    TGeoBBox *hbox = new TGeoBBox(.25, .1, .1);
    TGeoVolume *hbar = new TGeoVolume("hbar", hbox);
    solid->AddNode(hbar, 1, new TGeoTranslation(-.05, .4, 0.));
    solid->AddNode(hbar, 2, new TGeoTranslation(-.05, 0., 0.));

  } else if ( letter == "G" ) {
    // Add two vertical bar to the side of the letter
    TGeoBBox *vbox = new TGeoBBox(.1, .15, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.4, 0., 0.));

    // Add the two round parts to the top and bottom
    TGeoTubeSeg *tub = new TGeoTubeSeg(.15, .35, .1, 0., 180.);
    TGeoVolume *round = new TGeoVolume("round", tub);
    solid->AddNode(round, 1, new TGeoTranslation(-.15, .15, 0.));
    TGeoRotation *rot = new TGeoRotation("rot", 90., 0., 90., 270., 0., 0.);
    solid->AddNode(round, 2, new TGeoCombiTrans(-.15, -.15, 0., rot));

    // Add a stub to middle of the letter
    TGeoBBox *sbox = new TGeoBBox(.2, .1, .1);
    TGeoVolume *sbar = new TGeoVolume("sbar", sbox);
    solid->AddNode(sbar, 1, new TGeoTranslation(0., -.05, 0.));

  } else if ( letter == "H" ) {
    // Add the two vertical bars of the letter
    TGeoBBox *vbox = new TGeoBBox(.1, .5, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.4, 0., 0.));
    solid->AddNode(vbar, 1, new TGeoTranslation(.1, 0., 0.));

    // Add the horizontal bar of the letter
    TGeoBBox *hbox = new TGeoBBox(.15, .1, .1);
    TGeoVolume *hbar = new TGeoVolume("hbar", hbox);
    solid->AddNode(hbar, 1, new TGeoTranslation(-.15, 0., 0.));

  } else if ( letter == "I" ) {
    // Add the vertical bar of the letter
    TGeoBBox *vbox = new TGeoBBox(.1, .5, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.15, 0., 0.));

  } else if ( letter == "J" ) {
    // Add the vertical bar to the right of the letter
    TGeoBBox *vbox = new TGeoBBox(.1, .325, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(.1, 0.175, 0.));

    // Add the round part at the bottom
    TGeoTubeSeg *tub = new TGeoTubeSeg(.15, .35, .1, 180., 360.);
    TGeoVolume *round = new TGeoVolume("round", tub);
    solid->AddNode(round, 1, new TGeoTranslation(-.15, -.15, 0.));

  } else if ( letter == "K" ) {
    // Add the vertical bar to the left of the letter
    TGeoBBox *vbox = new TGeoBBox(.1, .5, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.4, 0., 0.));

    // Add the two diagonal of the letter
    TGeoArb8 *dbox = new TGeoArb8(.1);
    dbox->SetVertex(0, -.3, -.1);
    dbox->SetVertex(1, -.3, .1);
    dbox->SetVertex(2, .2, .5);
    dbox->SetVertex(3, .2, .3);
    dbox->SetVertex(4, -.3, -.1);
    dbox->SetVertex(5, -.3, .1);
    dbox->SetVertex(6, .2, .5);
    dbox->SetVertex(7, .2, .3);
    TGeoVolume *dbar = new TGeoVolume("diag", dbox);
    solid->AddNode(dbar, 1);
    TGeoRotation *rot = new TGeoRotation();
    rot->RotateX(180.);
    solid->AddNode(dbar, 2, rot);

  } else if ( letter == "L" ) {
    // Add the vertical bar of the letter
    TGeoBBox *vbox = new TGeoBBox(.1, .5, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.4, 0., 0.));

    // Add the horizontal bar of the letter
    TGeoBBox *hbox = new TGeoBBox(.25, .1, .1);
    TGeoVolume *hbar = new TGeoVolume("hbar", hbox);
    solid->AddNode(hbar, 1, new TGeoTranslation(-.05, -.4, 0.));

  } else if ( letter == "M" ) {
    // Add the vertical bar of the letter
    TGeoBBox *vbox = new TGeoBBox(.075, .5, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.425, 0., 0.));
    solid->AddNode(vbar, 2, new TGeoTranslation(.125, 0., 0.));

    // Add the diagonal bars
    TGeoArb8 *dbox = new TGeoArb8(.1);
    dbox->SetVertex(0, -.35, .3);
    dbox->SetVertex(1, -.35, .5);
    dbox->SetVertex(2, -.15, .1);
    dbox->SetVertex(3, -.15, -.1);
    dbox->SetVertex(4, -.35, .3);
    dbox->SetVertex(5, -.35, .5);
    dbox->SetVertex(6, -.15, .1);
    dbox->SetVertex(7, -.15, -.1);
    TGeoVolume *dbar = new TGeoVolume("diag", dbox);
    solid->AddNode(dbar, 1);
    TGeoRotation *rot = new TGeoRotation();
    rot->RotateY(180.);
    solid->AddNode(dbar, 2, new TGeoCombiTrans(-.3, 0., 0., rot));

  } else if ( letter == "N" ) {
    // Add the vertical bar of the letter
    TGeoBBox *vbox = new TGeoBBox(.075, .5, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.425, 0., 0.));
    solid->AddNode(vbar, 2, new TGeoTranslation(.125, 0., 0.));

    // Add the diagonal bar
    TGeoArb8 *dbox = new TGeoArb8(.1);
    dbox->SetVertex(0, -.35, .3);
    dbox->SetVertex(1, -.35, .5);
    dbox->SetVertex(2, .05, -.3);
    dbox->SetVertex(3, .05, -.5);
    dbox->SetVertex(4, -.35, .3);
    dbox->SetVertex(5, -.35, .5);
    dbox->SetVertex(6, .05, -.3);
    dbox->SetVertex(7, .05, -.5);
    TGeoVolume *dbar = new TGeoVolume("diag", dbox);
    solid->AddNode(dbar, 1);

  } else if ( letter == "O" ) {
    // Add two vertical bars to the sides of the letter
    TGeoBBox *vbox = new TGeoBBox(.1, .15, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.4, 0., 0.));
    solid->AddNode(vbar, 2, new TGeoTranslation(.1, 0., 0.));

    // Add the two round parts to the top and bottom
    TGeoTubeSeg *tub = new TGeoTubeSeg(.15, .35, .1, 0., 180.);
    TGeoVolume *round = new TGeoVolume("round", tub);
    solid->AddNode(round, 1, new TGeoTranslation(-.15, .15, 0.));
    TGeoRotation *rot = new TGeoRotation("rot", 90., 0., 90., 270., 0., 0.);
    solid->AddNode(round, 2, new TGeoCombiTrans(-.15, -.15, 0., rot));

  } else if ( letter == "P" ) {
    // Add the vertical bar of the letter
    TGeoBBox *vbox = new TGeoBBox(.1, .5, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.4, 0, 0.));

    // Add the two stubs to the right of the vertical bar
    TGeoBBox *sbox = new TGeoBBox(.1, .1, .1);
    TGeoVolume *sbar = new TGeoVolume("sbar", sbox);
    solid->AddNode(sbar, 1, new TGeoTranslation(-.2, 0., 0.));
    solid->AddNode(sbar, 2, new TGeoTranslation(-.2, .4, 0.));

    // Add the round part
    TGeoTubeSeg *tub = new TGeoTubeSeg(.1, .3, .1, -90., 90.);
    TGeoVolume *round = new TGeoVolume("round", tub);
    solid->AddNode(round, 1, new TGeoTranslation(-.1, .2, 0.));

  } else if ( letter == "Q" ) {
    // Add two vertical bars to the sides of the letter
    TGeoBBox *vbox = new TGeoBBox(.1, .15, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.4, 0., 0.));
    solid->AddNode(vbar, 2, new TGeoTranslation(.1, 0., 0.));

    // Add the two round parts to the top and bottom
    TGeoTubeSeg *tub = new TGeoTubeSeg(.15, .35, .1, 0., 180.);
    TGeoVolume *round = new TGeoVolume("round", tub);
    solid->AddNode(round, 1, new TGeoTranslation(-.15, .15, 0.));
    TGeoRotation *rot = new TGeoRotation("rot", 90., 0., 90., 270., 0., 0.);
    solid->AddNode(round, 2, new TGeoCombiTrans(-.15, -.15, 0., rot));

    // Add the diagonal bar
    TGeoArb8 *dbox = new TGeoArb8(.1);
    dbox->SetVertex(0, -.15, -.2);
    dbox->SetVertex(1, -.15, 0.);
    dbox->SetVertex(2, .2, -.3);
    dbox->SetVertex(3, .2, -.5);
    dbox->SetVertex(4, -.15, -.2);
    dbox->SetVertex(5, -.15, 0.);
    dbox->SetVertex(6, .2, -.3);
    dbox->SetVertex(7, .2, -.5);
    TGeoVolume *dbar = new TGeoVolume("diag", dbox);
    solid->AddNode(dbar, 1);

  } else if ( letter == "R" ) {
    // Add the vertical bar of the letter
    TGeoBBox *vbox = new TGeoBBox(.1, .5, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.4, 0, 0.));

    // Add the two stubs to the right of the vertical bar
    TGeoBBox *sbox = new TGeoBBox(.1, .1, .1);
    TGeoVolume *sbar = new TGeoVolume("sbar", sbox);
    solid->AddNode(sbar, 1, new TGeoTranslation(-.2, 0., 0.));
    solid->AddNode(sbar, 2, new TGeoTranslation(-.2, .4, 0.));

    // Add the round part
    TGeoTubeSeg *tub = new TGeoTubeSeg(.1, .3, .1, -90., 90.);
    TGeoVolume *round = new TGeoVolume("round", tub);
    solid->AddNode(round, 1, new TGeoTranslation(-.1, .2, 0.));

    // Add the diagonal bar
    TGeoArb8 *dbox = new TGeoArb8(.1);
    dbox->SetVertex(0, -.3, -.1);
    dbox->SetVertex(1, -.1, -.1);
    dbox->SetVertex(2, .2, -.5);
    dbox->SetVertex(3, 0., -.5);
    dbox->SetVertex(4, -.3, -.1);
    dbox->SetVertex(5, -.1, -.1);
    dbox->SetVertex(6, .2, -.5);
    dbox->SetVertex(7, 0., -.5);
    TGeoVolume *dbar = new TGeoVolume("diag", dbox);
    solid->AddNode(dbar, 1);

  } else if ( letter == "S" ) {
    // Add the two round parts to the top left and bottom right
    TGeoTubeSeg *tub = new TGeoTubeSeg(.1, .3, .1, 90., 270.);
    TGeoVolume *round = new TGeoVolume("round", tub);
    solid->AddNode(round, 1, new TGeoTranslation(-.2, .2, 0.));
    TGeoRotation *rot = new TGeoRotation();
    rot->RotateZ(180);
    solid->AddNode(round, 2, new TGeoCombiTrans(-.1, -.2, 0., rot));

    // Add the two round parts to the top right and bottom left
    TGeoTubeSeg *tub2 = new TGeoTubeSeg(.1, .3, .1, 0., 90.);
    TGeoVolume *round2 = new TGeoVolume("round2", tub2);
    solid->AddNode(round2, 1, new TGeoTranslation(-.1, .2, 0.));
    solid->AddNode(round2, 2, new TGeoCombiTrans(-.2, -.2, 0., rot));

    // Add wedges to connect the two round parts
    TGeoBBox *wbox = new TGeoBBox(.05, .1, .1);
    TGeoVolume *wbar = new TGeoVolume("wbar", wbox);
    solid->AddNode(wbar, 1, new TGeoTranslation(-.15, 0., 0.));
    solid->AddNode(wbar, 1, new TGeoTranslation(-.15, .4, 0.));
    solid->AddNode(wbar, 1, new TGeoTranslation(-.15, -.4, 0.));

  } else if ( letter == "T" ) {
    // Add the vertical bar of the letter
    TGeoBBox *vbox = new TGeoBBox(.1, .4, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.15, -.1, 0.));

    // Add the horizontal bar of the letter
    TGeoBBox *hbox = new TGeoBBox(.35, .1, .1);
    TGeoVolume *hbar = new TGeoVolume("hbar", hbox);
    solid->AddNode(hbar, 1, new TGeoTranslation(-.15, .4, 0.));

  } else if ( letter == "U" ) {
    // Add the vertical bars of the letter
    TGeoBBox *vbox = new TGeoBBox(.1, .325, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.4, 0.175, 0.));
    solid->AddNode(vbar, 2, new TGeoTranslation(.1, 0.175, 0.));

    // Add the round part at the bottom
    TGeoTubeSeg *tub = new TGeoTubeSeg(.15, .35, .1, 180., 360.);
    TGeoVolume *round = new TGeoVolume("round", tub);
    solid->AddNode(round, 1, new TGeoTranslation(-.15, -.15, 0.));

  } else if ( letter == "V" ) {
    // Add the two diagonals of the letter
    TGeoArb8 *dbox = new TGeoArb8(.1);
    dbox->SetVertex(0, -.3, .5);
    dbox->SetVertex(1, -.5, .5);
    dbox->SetVertex(2, -.25, -.5);
    dbox->SetVertex(3, -.05, -.5);
    dbox->SetVertex(4, -.3, .5);
    dbox->SetVertex(5, -.5, .5);
    dbox->SetVertex(6, -.25, -.5);
    dbox->SetVertex(7, -.05, -.5);
    TGeoVolume *dbar = new TGeoVolume("diag", dbox);
    solid->AddNode(dbar, 1);
    TGeoRotation *rot = new TGeoRotation();
    rot->RotateY(180.);
    solid->AddNode(dbar, 2, new TGeoCombiTrans(-.3, 0., 0., rot));

  } else if ( letter == "W" ) {
    // Add the two outside bars of the letter
    TGeoArb8 *dbox = new TGeoArb8(.1);
    dbox->SetVertex(0, -.35, .5);
    dbox->SetVertex(1, -.5, .5);
    dbox->SetVertex(2, -.4, -.5);
    dbox->SetVertex(3, -.25, -.5);
    dbox->SetVertex(4, -.35, .5);
    dbox->SetVertex(5, -.5, .5);
    dbox->SetVertex(6, -.4, -.5);
    dbox->SetVertex(7, -.25, -.5);
    TGeoVolume *dbar = new TGeoVolume("diag", dbox);
    solid->AddNode(dbar, 1);
    TGeoRotation *rot = new TGeoRotation();
    rot->RotateY(180.);
    solid->AddNode(dbar, 2, new TGeoCombiTrans(-.3, 0., 0., rot));

    // Add the two inside bars of the letter (slighty different tilt angle)
    TGeoArb8 *dboxc = new TGeoArb8(.1);
    dboxc->SetVertex(0, -.25, -.5);
    dboxc->SetVertex(1, -.4, -.5);
    dboxc->SetVertex(2, -.225, .5);
    dboxc->SetVertex(3, -.075, .5);
    dboxc->SetVertex(4, -.25, -.5);
    dboxc->SetVertex(5, -.4, -.5);
    dboxc->SetVertex(6, -.225, .5);
    dboxc->SetVertex(7, -.075, .5);
    TGeoVolume *dbarc = new TGeoVolume("diagc", dboxc);
    solid->AddNode(dbarc, 1);
    solid->AddNode(dbarc, 2, new TGeoCombiTrans(-.3, 0., 0., rot));

  } else if ( letter == "X" ) {
    // Add the two diagonals of the letter
    TGeoArb8 *dbox = new TGeoArb8(.1);
    dbox->SetVertex(0, -.3, -.5);
    dbox->SetVertex(1, -.5, -.5);
    dbox->SetVertex(2, 0., .5);
    dbox->SetVertex(3, .2, .5);
    dbox->SetVertex(4, -.3, -.5);
    dbox->SetVertex(5, -.5, -.5);
    dbox->SetVertex(6, 0., .5);
    dbox->SetVertex(7, .2, .5);
    TGeoVolume *dbar = new TGeoVolume("diag", dbox);
    solid->AddNode(dbar, 1);
    TGeoRotation *rot = new TGeoRotation();
    rot->RotateY(180.);
    solid->AddNode(dbar, 2, new TGeoCombiTrans(-.3, 0., 0., rot));

  } else if ( letter == "Y" ) {
    // Add the two diagonals of the letter first
    TGeoArb8 *dbox = new TGeoArb8(.1);
    dbox->SetVertex(0, -.3, .5);
    dbox->SetVertex(1, -.5, .5);
    dbox->SetVertex(2, -.25, 0.);
    dbox->SetVertex(3, -0.05, 0.);
    dbox->SetVertex(4, -.3, .5);
    dbox->SetVertex(5, -.5, .5);
    dbox->SetVertex(6, -.25, 0.);
    dbox->SetVertex(7, -0.05, 0.);
    TGeoVolume *dbar = new TGeoVolume("diag", dbox);
    solid->AddNode(dbar, 1);
    TGeoRotation *rot = new TGeoRotation();
    rot->RotateY(180.);
    solid->AddNode(dbar, 2, new TGeoCombiTrans(-.3, 0., 0., rot));

    // Add the vertical bar of the letter
    TGeoBBox *vbox = new TGeoBBox(.1, .25, .1);
    TGeoVolume *vbar = new TGeoVolume("vbar", vbox);
    solid->AddNode(vbar, 1, new TGeoTranslation(-.15, -.25, 0.));

  } else if ( letter == "Z" ) {
    // Add the horizontal bars of the letter
    TGeoBBox *hbox = new TGeoBBox(.35, .1, .1);
    TGeoVolume *hbar = new TGeoVolume("hbar", hbox);
    solid->AddNode(hbar, 1, new TGeoTranslation(-.15, .4, 0.));
    solid->AddNode(hbar, 2, new TGeoTranslation(-.15, -.4, 0.));

    // Add the diagonal bar
    TGeoArb8 *dbox = new TGeoArb8(.1);
    dbox->SetVertex(0, -.3, -.3);
    dbox->SetVertex(1, -.5, -.3);
    dbox->SetVertex(2, 0., .3);
    dbox->SetVertex(3, .2, .3);
    dbox->SetVertex(4, -.3, -.3);
    dbox->SetVertex(5, -.5, -.3);
    dbox->SetVertex(6, 0., .3);
    dbox->SetVertex(7, .2, .3);
    TGeoVolume *dbar = new TGeoVolume("diag", dbox);
    solid->AddNode(dbar, 1);

  } else if ( letter == " " ) {
    // If the letter is a blank space, don't fill it with anything

  } else {
    std::cerr << "Letter " << letter << " currently not supported" << std::endl;
  }
  return solid;
}

int main(int argc, char ** argv) {

  // Read the input string to be produced from the command line, make it upper case
  std::string string = argv[1];
  std::transform(string.begin(), string.end(), string.begin(),
                 [](unsigned char c) { return std::toupper(c); });
  size_t nchars = string.size();

  // Initialize the geometry manager and the parent volume according to the string length
  double dx(nchars*.5), dy(.5), dz(.1);
  TGeoManager *geom = new TGeoManager("top", "Top volume");
  TGeoBBox *topbox = new TGeoBBox(dx, dy, dz);
  TGeoVolume *top = new TGeoVolume("topbox", topbox);
  top->SetVisibility(kFALSE);
  geom->SetTopVolume(top);

  // Add letters to the parent volume one by one
  std::string buff;
  double posx;
  size_t i;
  for (i = 0; i < nchars; i++) {
    buff = string[i];
    TGeoVolume* vol = Letter(buff);
    posx = (1.-nchars)*.5 + i;
    top->AddNode(vol, 1, new TGeoTranslation(posx, 0., 0.));
  }

  // To create a sample of points inside the letter(s), one needs to generate random points
  // inside of the bounding box and only keep the ones inside the letter(s) volume(s)
  TRandom3 rdmzer(time(NULL));
  size_t nsamples = 1000;
  if ( argc > 2 )
      nsamples = atoi(argv[2]);
  std::vector<std::vector<double>> points;

  std::vector<double> temp;
  TGeoNode *node(NULL);
  top->VisibleDaughters(kTRUE);
  top->Draw();
  i = 0;
  while (i < nsamples) {
    // First generate a point randomly insine the box
    temp = {-dx+2*dx*rdmzer.Rndm(), -dy+2*dy*rdmzer.Rndm(), -dz+2*dz*rdmzer.Rndm()};
    geom->SetCurrentPoint(&temp[0]);

    // Check if it is inside any of the letter nodes
    node = geom->FindNode();
    if ( !node )
	continue;
    if ( !node->IsOnScreen() )
	continue;

    // If it is inside, add it to the samples
    points.push_back(temp);
    i++;
  }

  // Write the points to a TFile for future use
  std::string filename = "samples_"+string+".root";
  TFile f(filename.c_str(), "RECREATE");
  TTree *tpoints = new TTree;
  std::vector<double> *point(NULL);
  tpoints->Branch("points", &point);
  for (i = 0; i < points.size(); i++) {
    point = new std::vector<double>(points[i]);
    tpoints->GetBranch("points")->Fill();
    delete point;
  }
  tpoints->Write("tree");
  f.Close();

  // Draw the letters
  TCanvas *c = new TCanvas("c", "c", 1600, 800);
//  top->RandomPoints(nsamples);
  top->Draw("");
  gPad->GetView()->TopView();
  c->SaveAs(std::string("test_"+string+".pdf").c_str());
  c->SaveAs(std::string("test_"+string+".root").c_str());
}
