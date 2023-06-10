#ifndef NODE_H_
#define NODE_H_

#include <iostream>
#include <stdlib.h>
#include <string>
#include <string>
#include <fstream>

#include "matrices/matrix.h"      // for CP matrix

using namespace std;   // allows direct reference to objects such as cout

class node
{
protected://attributes
   int id;                         // node ID number
   int height;                     // height of node is ....
   int anc;                        // ancestor of the node
   int left;                       // left descendant of the node
   int right;                      // right descendant of the node
   int row_index;                  // index of row in data matrix
   int branch_class;               // Used by branch-class models; holds a "branch label"
   int clade_class;                // Used by non-homogeneous models; holds a "clade label"
   string state;                   // character state at the node
   string label;                   // label = tip name; otherwise it should be "internal"
   double blength;                 // length of branch to ancestor!!!!

public:// constructor
   node ();

public://attributes
   matrix  CP_mat;               // conditional prob matrix; reallocate to 4/20/61 x num_patterns; could move this to a model object in future too? (but keeping it in node keeps track of context!)
   matrix  TP_mat;               // Trying to abandon the derived type matrix in a node and only use the base class matrix; we shall see how it turns out

   // ***********//
   // refactor so that DNA_matrix is unnecessary (I think this is done!)
   //DNA_matrix* TP_matrix;        // ptr to a transition prob matrix; type and dim will depend on AA, DNA, or CODON data; lives in a model object?
   //DNA_matrix DNA_mat;           //

public:// functions
   void set_id(int node_num);              // set ID number of a node object
   void set_anc(int node_num);             // set ANCESTOR of a node object
   void set_left(int node_num);            // set LEFT-DESCENDANT of a node object
   void set_right(int node_num);           // set RIGHT-DESCENDANT of a node object
   void set_row_index(int row_num);        // set the row index for TIP NODES ONLY
   void set_state(string input_state);     // set the character state at a node
   void set_label(string input_str);       // set the node label
   void set_blength(double input_length);  // set branch length to ancestor of present node
   void set_height(int ht);                // set the height value of a node
   //void set_TPmatrix(DNA_matrix* ptr);     // set ptr to the a derived matrix object
   void set_branch_class(int id);          // set the node id for a particular class of a "branch class modes"
   void set_clade_class(int id);           // set the node id for a particular clade class

   //void init_DNA_CPmat(int pattern);       // Moved to Likelihood engine why should node have to keep track of something to do with calcualting the likelihood given a data type!!!

   int get_id();                           // returns ID of the current node
   int get_anc();                          // returns ID of ancestor of present node
   int get_left();                         // returns ID of left descendant of present node
   int get_right();                        // returns ID of right descendant of present node
   int get_row_index();                    // returns index of row in data matrix for the involved OTU
   int get_height();                       // returns the height of a node
   int get_branch_class();                 // returns the node id for a particular class of a "branch class modes"
   int get_clade_class();                  // returns the node id for a particular clade-class
   string get_label();                     // returns the node label
   string get_state();                     // returns the state of the node
   bool is_tip();                          // returns 1 if the node is a tip; else 0
   bool is_root();                         // returns 1 if the node is the root; else 0
   void show();                            // prints the node attributes to the screen
   double get_blength();                   // returns length of branch from present node to its ancestor
};

#endif /*NODE_H_*/

