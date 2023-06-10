#include <iostream>
#include <stdlib.h>
#include <string>
#include <algorithm>
#include <fstream>
#include "node.h"

using namespace std;   // allows direct reference to objects such as cout

//------------------------------------------------------------------------------
// CONSTRUCTOR:
//------------------------------------------------------------------------------
node::node():id(-1), height(-1), anc(-1), left(-1), right(-1), row_index(-1),
branch_class(0), clade_class(0), state("NULL"), label("NULL"), blength(-1)
{
   //set_state("null");
   //set_label("null");
}

//------------------------------------------------------------------------------
// public: show all data for a node
//------------------------------------------------------------------------------
void node::show()
{
   cout << "label             : " << label << "\n";
   cout << "node ID           : " << id << "\n";
   cout << "branch_class      : " << branch_class << "\n";
   cout << "branch length     : " << blength << "\n";
   cout << "row index         : " << row_index << "\n";
   cout << "state             : " << state << "\n";
   cout << "left descendant   : " << left << "\n";
   cout << "right descendant  : " << right << "\n";
   cout << "ancestor          : " << anc << "\n";
}

//------------------------------------------------------------------------------
// public: set ID number
//------------------------------------------------------------------------------
void node::set_id(int node_num)
{
   id = node_num;
}

//------------------------------------------------------------------------------
// public: set ancestor of node
//------------------------------------------------------------------------------
void node::set_anc(int node_num)
{
   anc = node_num;
}

//------------------------------------------------------------------------------
// public: set left descendant of node
//------------------------------------------------------------------------------
void node::set_left(int node_num)
{
   left = node_num;
}

//------------------------------------------------------------------------------
// public: set right descendant of node
//------------------------------------------------------------------------------
void node::set_right(int node_num)
{
   right = node_num;
}

//------------------------------------------------------------------------------
// public: set number of row of involved OTU in the se data matrix
//------------------------------------------------------------------------------
void node::set_row_index(int row_num)
{
   row_index = row_num;
}

//------------------------------------------------------------------------------
// public: set the character state of the node; default value = NULL
//------------------------------------------------------------------------------
void node::set_state(string input_state)
{
   state = input_state;
}

//------------------------------------------------------------------------------
// public: set the node label; default value = NULL; interior nodes remain = NULL
//------------------------------------------------------------------------------
void node::set_label(string input_str)
{
   label = input_str;
}

//------------------------------------------------------------------------------
// public: set the node label; default value = NULL; interior nodes remain = NULL
//------------------------------------------------------------------------------
void node::set_blength(double input_length)
{
   blength = input_length;
}

//------------------------------------------------------------------------------
// public: set the node label; default value = -1
//------------------------------------------------------------------------------
void node::set_height(int ht)
{
	height = ht;
}

//------------------------------------------------------------------------------
// public:  set the node id for a particular class of a "branch class modes"
//------------------------------------------------------------------------------
void node::set_branch_class(int id)
{
   branch_class = id;
}


//------------------------------------------------------------------------------
// public:  set the node id for a particular clade class
//------------------------------------------------------------------------------
void node::set_clade_class(int id)
{
   clade_class = id;
}

//------------------------------------------------------------------------------
// public: get ID
//------------------------------------------------------------------------------
int node::get_id()
{
   return id;
}

//------------------------------------------------------------------------------
// public: get ancestral node
//------------------------------------------------------------------------------
int node::get_anc()
{
   return anc;
}

//------------------------------------------------------------------------------
// public: get left descendant
//------------------------------------------------------------------------------
int node::get_left()
{
   return left;
}

//------------------------------------------------------------------------------
// public: get right descendant
//------------------------------------------------------------------------------
int node::get_right()
{
   return right;
}

//------------------------------------------------------------------------------
// public: get branch length
//------------------------------------------------------------------------------
double node::get_blength()
{
   return blength;
}

//------------------------------------------------------------------------------
// public: get node label
//------------------------------------------------------------------------------
string node::get_label()
{
   return label;
}

//------------------------------------------------------------------------------
// public: get node state
//------------------------------------------------------------------------------
string node::get_state()
{
   return state;
}

//------------------------------------------------------------------------------
// public: get index of row in data matrix for the involved OTU
//------------------------------------------------------------------------------
int node::get_row_index()
{
   return row_index;
}

//------------------------------------------------------------------------------
// public: get height of a node
//------------------------------------------------------------------------------
int node::get_height()
{
	return height;
}

//------------------------------------------------------------------------------
// public: returns the node id for a particular class of a "branch class modes"
//------------------------------------------------------------------------------
int node::get_branch_class()
{
   return branch_class;
}


//------------------------------------------------------------------------------
// public: returns the node id for a particular clade-class
//------------------------------------------------------------------------------
int node::get_clade_class()
{
   return clade_class;
}


//------------------------------------------------------------------------------
// public:
//------------------------------------------------------------------------------
bool node::is_tip()
{
   if(left == -1 && right == -1) {return true;}
   else {return false;}
}


//------------------------------------------------------------------------------
// public:
//------------------------------------------------------------------------------
bool node::is_root()
{
   // think about this later; root node should always have ID==0 as well!
   if(anc == -1) {return true;}
   else {return false;}
}




