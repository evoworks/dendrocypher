#ifndef FORESTRY_H_
#define FORESTRY_H_

#include <iostream>
#include <string>
#include <vector>
#include "tree.h"

using namespace std;

class forestry
{
private://attributes
	int num_trees;                     // number of trees in the tree file
	string tree_file;                  // name of tree file; read from options.txt
public:  //TEMP PUBLIC
	tree* my_forest;                   // base tree obj that serves to hold a vector of all the tree objects called "population"
	vector<string> newick_trees;       // each string in vector is a newick tree

private://methods
	bool validate(string newick);      // called by "read_newick_trees"; check format of the newick string: (i) matched parens, (2) ....

public://constructors
	forestry ();                       // default: do nothing
	forestry (string file, tree* ptr); // constructor is defined in forestry.cpp; file=treefile

public://methods
	int get_num_trees();               // returns number of trees
	void set_forest(tree* ptr);        // sets an attribute that points to a tree object that contains *the population of trees*
	void set_treefile(string file);    // use this to set the name of the file holding the newick trees
	void read_newick_trees();          // read newick trees and put text in an private vector called "newick_trees"
	void show_newick_trees();          // print the contents of the private newick_trees array
	string get_newick_tree(int id);    // return newick tree at element i in the "newick_trees" array (this is NOT from the BTM)
	void populate_forest();            // create new tree objects in a static vector in the base tree class that holds a "population" of trees
	void create_tree_objects();        // loop over newick tree array array and create a tree object for each newick tree
	void show_tree(int tree_id);       // print tree to screen as it exists as a BTM!!!
	void check_trees();                // to check BTM, print newick string and print euler tour; should be same
	void delete_trees();               // see note in the implementation: start at tree i=1!
	string get_newick(int tree_id);    // return a newick string for the current BTM and its blengths (if any)
	tree* get_tree(int i);             // return a pointer to tree_i; UNTESTED in forestry!!!
};

#endif /* FORESTRY_H_ */

