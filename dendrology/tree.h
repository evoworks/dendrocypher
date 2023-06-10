#ifndef TREE_H_
#define TREE_H_

#include <iostream>
#include <sstream>
#include <stdlib.h>
#include <string>
#include <fstream>
#include <vector>
#include "node.h"

using namespace std;

class tree
{
protected://attributes
    string newick_tree;              // newick tree string
    stringstream mySS;

    // fix this max OTU limitation!!!
    char otu[100][100];              // limits on OUTs: max 30 OTUs, and max 100 char label (re-factor using a vector? YUK.)
    int num_otu;                     // number of OTUs
    int num_nodes;                   // number of nodes in the tree
    int num_branches;                // number of branches in the tree
    int the_root;                    // ID number of the root node; should = 0
    int tree_height;                 // Height of a tree = height of the root node
    int sum_stats[4];                // OUTDATED !?!?!?!?!
    int index;                       // index =???; LOOK this up!!!; perhaps us an ID, and index to location in the population vector
    bool rooted;                     // true if basal DICHOTOMY; unrooted achieves TRICHOTOMY at root node only by FIXING branch length to left_of_root = 0; root node has same id!

public: //attributes
    static vector<tree*> population;     // defined as a static variable in the implementation file; static ensures that one copy is shared by all tree objects

protected://methods
    void uppass(node &input_node);                       // downpass recursion on left & right descendants
    void euler(int i);                                   // algorithm to get parenthetical notation (prints Newick tree to screen via cout)
    void labeled_euler(int i, int flag);                 // shows Newick tree with BRANCH CLASS ID values
    string euler2(int i, string str);                    // returns an string with tree in parenthetical notation
    int get_right_descendant(int q, string local_tree);  // imp fix to the creat_BTM function; add more notes later
    void set_node_heights();                             // used to call the recursive algorithm below; sets height of each node from root
    void recursive_heights(int ht, node &input_node);    // recursive algorithm that sets the height values for all ancestors of the input node

public://constructors
    tree();
    tree(string input_str);
    ~tree();

private://attributes
    vector<int> traversal_order;           // traversal order of nodes; downpass, start at node 0; for uppass, start at num-nodes-1; init w/ -1 in "populate_nodes()"
public://attributes
    vector<node> nodes;                    // a populate method puts the req'd num of generic nodes in the vector; create_BTM2 does the real work

public://methods
    int get_num_nodes();                    // returns number of nodes in the tree
    int get_num_otu();                      // returns number of OTUs (tip nodes) in the tree
    int get_downpass_node(int rank);        // returns the node ID for a downpass tree traversal given the current place (i.e., rank) in the traversal process (0, 1, 2, ... num nodes)
    void stats();                           // print some basic summary tree stats to the screen
    void set_newick(string input_str);      // the variable "newick_tree" is set to a string = newick formatted tree
    void set_num_nodes();                   // num_otu = 2 x number of commas + 1 in the newick string; req'd for populate_nodes() below.
    void populate_nodes();                  // place a node object in vector for each node in the tree
    void create_btm();                      // construct a binary tree model(BTM)from a newick tree
    void create_btm2();                     // this version fixes a bug in the above method!!!!
    void load_blengths();                   // does what it says on the tin
    int load_branch_marks();                // load branch marks identified within a newick tree by the "#" symbol
    int load_clade_marks();                 // load clade marks identified within a newick tree by the "$" symbol
    void show_nodes();                      // show nodes of tree; too much for a big tree!!!
    void traverse();                        // this function calls private function "downpass"
    void print_downpass();                  // print order of nodes in a DOWNPASS to the screen
    void print_uppass();                    // print order of nodes in a UPPASS to the screen
    void show_heights();                    //  ...ADD SOME NOTES ...
    void show_branch_classes();             // show to screen a variety of info related to branch classes

    string get_node_label(int node_num);    // supply node ID and the node label is returned
    int get_root();                         // should ALWAYS be zero
    void euler_tour();                      // initiate the euler tour at the root: overloaded version below returns string
    void labeled_euler_tour(int flag);      // shows newick tree to screen with BRANCH CLASS ID values
    string euler_tour(string str);          // now doing nothing with str; in future could use to hold a filename
    void delete_trees();                    // delete tree objects from heap to prevent memory leak
    void show_states(bool internal);        // show node states; if internal = 0, don't show; if internal = 1, show states
    void deroot();                          // convert a rooted tree to an unrooted tree by fixing the blength to left_of_root = 0; root node ID stays same; bool rooted is set to false
    bool is_rooted();                       // returns true or false depending on status or "rooted" variable
    void set_root(int target, bool blengths, bool bmarks);  // ...ADD SOME NOTES ...
    bool compare_trees_from(int id);        // compare trees in population starting from tree[0] to tree[id]; check if identical and return true (false if topologies differ); does NOT compare b.l.s or marks!
    // void set_root(string label);
    // set_basal_trichotomy                 // Yields unrooted tree; should code so that it helps to deal with branch lengths during the optimization!!!!

    
public://treeCypher methods
    void show_numbered_tree();              // show tree with ALL nodes numbered
    void show_marked_tree();                // show tree with any branch marks, or clade marks, if present in tree data structure
    void show_branch_class_table();         // show table of branch classes on screen (mostly for de-bugging)
    void print_numbered_tree();             // print tree to a file with ALL nodes numbered    void print_bipartitions(int id);        // print birpartions to a file for tree at index = id
    void print_marked_tree();               // print tree AND include any branch or clade marks
    
    void show_node_marks(int node);                // shows any marks (branch or clade) or labels for the node.
    void set_branch_mark(int node, int mark);      // add the mark (integer) to the specified node
    void remove_branch_mark(int node);             // remove the mark (integer) to the specified node
    void numbered_euler(int i, bool print);        // maybe make this private?
    void marked_euler(int i, bool print);          // maybe make this private?
    void list_descendants(bool show, bool print);  //
    void visit_descendants(int i);                 //

public://friends
    friend class forestry;                // forestry obj can access tree private parts
};

#endif /*TREE_H_*/


