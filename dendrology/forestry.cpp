#include <iostream>
#include <stdlib.h>
#include <string>
#include <sstream>
#include "forestry.h"

using namespace std;

//------------------------------------------------------------------------------
// FORESTRY: constructor
//------------------------------------------------------------------------------
forestry::forestry():num_trees(-1), tree_file("")
{
   //newick_trees.reserve(100);
}


//------------------------------------------------------------------------------
// FORESTRY: constructor
//------------------------------------------------------------------------------
forestry::forestry(string file, tree* ptr):num_trees(-1)
{
   tree_file = file;            // name of file that contains a list of newick formatted trees
   my_forest = ptr;             // set ptr to a base tree object that serves to hold the collection of tree objs in a vector called "population"
   read_newick_trees();         // read the newick tree strings from the file and "push_back" onto the vector "newick_trees"
   populate_forest();           // populates a static vector of trees contained by a base class tree object pointed to by the "my_forest" ptr; see notes in the implementation file; tree constructor does the real work!
   create_tree_objects();       // Tree objects converted to actual BTM based on the newick string in the "newick_trees" vector
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
int forestry::get_num_trees()
{
   return num_trees;
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void forestry::set_forest(tree* ptr)
{
	my_forest = ptr;
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void forestry::read_newick_trees()
{
   int index=0;
   string input_line;

   ifstream file_in(tree_file.c_str());   // c_str() method allows C-style string arg for ifstream!
   if(!file_in)
   {
      cout << "The file " << tree_file << " could NOT be opened.\n\n";
      num_trees=0;
      cout << "Sorry, this application must terminate.\n\n";
      system("PAUSE");
      terminate();
   }

   cout << "\tReading file: " << tree_file << "\n";
   while (!file_in.eof())
   {
         getline(file_in, input_line);
         if(!input_line.empty()) //only work on non-empty lines
         {
         	bool check=validate(input_line);
         	if(check == true)
         	{
         		//newick_trees.push_back(index) = input_line;
         		newick_trees.push_back(input_line);
         		index++;
         	}
         	else
         	{
         		cout << "\t\tFailed validation: skipping tree " << index+1 << "\n";
         	}
         }
   }
   file_in.close();  // close the file.
   num_trees=index;  // set the number of trees in this object
   //cout << "\n";
}


//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void forestry::show_newick_trees()
{
   cout << "\tNumber of trees = " << num_trees << "\n";
   cout << "\tCurrent trees:\n";
   cout << "\t--------------\n";
   for (int i=0; i < num_trees; i++)
   {
      cout << "\tTree " << (i+1) << " = " << newick_trees[i] << "\n";
   }
   cout << "\t--------------\n";
}

//------------------------------------------------------------------------------
// PUBLIC: return newick tree at element i in the "newick_trees" array (this
//         is NOT from the BTM)
//------------------------------------------------------------------------------
string forestry::get_newick_tree(int id)
{
	return newick_trees[id];
}


//------------------------------------------------------------------------------
// PUBLIC: create new tree objects in a static vector in the base tree class that
// holds a "population" of trees
//
// NOTE: The "my_forest" member variable is a ptr that points to a tree object!
//       This should point to the initial tree object (allocated in the registry, or
//       elsewhere if using the options file) that provides access to the tree population.
//       The population is defined as a static variable in the tree implementation file;
//       static ensures that one copy is shared by all tree objects. Thus forestry can
//       manage the allocation of new trees by pushing their pointers onto the back
//       of a single (shared) tree population vector at runtime.
//------------------------------------------------------------------------------
void forestry::populate_forest()
{
	// NOTE 1: Previously, constructing a base class tree in the "registry" to function as a
	// container for the vector of trees that get populated by this function yielded a pointer
	// to a single tree at location = 0. This means than an "extra" tree object will end up
	// in the "population" vector. To avoid this, the populate method stop(ed) as (num_trees - 1).

	// NOTE 2: (18 Feb 2014) To allow this method to work for BOTH the script interface and
	// the option file interface, I stopped the registry from loading a pointer to a tree in
	// the population vector.  This allows both interfaces to work with the code that is
	// currently implemented below. The important advantage is that the code no longer needs
	// to be re-compiled to switch between the two user interfaces!
    
    // NOTE 3: The base class tree object is construcetd in main.cpp in the program DendroCypher,
    // which also uses this code.


	//cout << "forestry::populate_forest: num_trees = " << num_trees << "\n";
	//cout << "forestry::populate_forest: tree popn size = " << my_forest->population.size() << "\n";

	//for (int i=0; i < (num_trees-1); i++)  /****** OLD: formerly used for script interface ******/
	for (int i=0; i < (num_trees); i++)
	{
		// Below, it pushes a ptr to the new tree obj (on the heap) into the tree
		// population vector. Thus, all that is needed is the key word "new"

		// OLD: my_forest->population.push_back(new tree);  /****** OLD: formerly used for script interface ******/

		if(i==0) my_forest->population.push_back(my_forest); // First element always gets pointer (my-forest) to the tree that was allocated in the registry (or elsewhere is using the options file; or main in DendroCypher)
		else my_forest->population.push_back(new tree);      // The rest of the trees must be newly allocated here
	}

	// NOTE: this method is separate from "create_tree_objects" below because I might want
	// to have a set of generic tree objects that get their structure from someplace other than a
	// file full of newick tree strings
	//
	//cout << "forestry::populate_forest: num_trees = " << num_trees << "\n";
	//cout << "forestry::populate_forest: tree popn size = " << my_forest->population.size() << "\n";
}


//------------------------------------------------------------------------------
// PUBLIC: use this to set the "tree_file"
//------------------------------------------------------------------------------
void forestry::set_treefile(string file)
{
   tree_file = file;
}


//------------------------------------------------------------------------------
// public:
//------------------------------------------------------------------------------
void forestry::create_tree_objects()
{
   cout << "\tNumber of trees = " << num_trees << "\n\n";
   //cout << "forestry::create_tree_objects: tree popn size = " << my_forest->population.size() << "\n\n";
   for (int i=0; i < num_trees; i++)
   {
      cout << "\tLoading tree " << i+1 << ": ";

      //cout << "forestry vector: " << newick_trees[i] << "\n";
      my_forest->population[i]->set_newick(newick_trees[i]);   // copy newick tree to a string variable in the tree obj
      //cout << "tree obj string: " << my_forest->population[i]->newick_tree << "\n";

      my_forest->population[i]->set_num_nodes();               // read newick string, then set number of nodes
      my_forest->population[i]->populate_nodes();              // populate a vector of nodes with node objects

      //int right = trees[i].get_right_descendant(8);
      //cout << "confirming that the right descendant = " << right << "\n";

      cout << " ...creating BTM";
      my_forest->population[i]->create_btm2();

       /*
      //For debugging...
      cout << "\n\n";
      my_forest->population[i]->show_nodes();
      //...For debugging
        */

      cout << " ...traversing";
      my_forest->population[i]->traverse();

      // For debugging...
      //cout << "\n\n";
      //my_forest->population[i]->print_downpass();
      //

      cout << " ...loading branch lengths";
      my_forest->population[i]->load_blengths();
      cout << " ...done.\n";

      my_forest->population[i]->set_node_heights();
      my_forest->population[i]->show_heights();

      // For debugging...
      //cout << "\n\n";
      //my_forest->population[i]->print_downpass();
      //

      //cout << "\n";

      cout << "\tNewick " << (i+1) << " = " << newick_trees[i] << "\n";
      cout << "\tBTM    " << (i+1) << " = ";
      my_forest->population[i]->euler_tour();
      cout << "\n\n";

      //my_forest->population[i]->stats();
      //my_forest->population[i]->show_nodes();
      //my_forest->population[i]->print_downpass();
      //my_forest->population[i]->print_uppass();
      //system("PAUSE");
      //cin.get();
   }
   //cout << "\n";
   //cout << "forestry::create_tree_objects: Number of trees = " << num_trees << "\n\n";
   //cout << "forestry::create_tree_objects: tree popn size = " << my_forest->population.size() << "\n\n";
}

//------------------------------------------------------------------------------
// private
//------------------------------------------------------------------------------
bool forestry::validate(string newick)
{
    //cout << "Calling the validate function on:\n";
    //cout << newick;

    //check for ";" at the end of the line
    int size = newick.size();                         // size counts from one
    int last = size-1;                                // position counts from zero; subtract 1 to get last pos
    int pos = newick.find_first_of(";");              // find pos of first semicolon (should be 1 from last)
    if( pos == (last-1) || pos == last)
    {
    	// do nothing
    }
    else
    {
    	cout << "\t\tdataset::validate::ERROR...semicolon?\n";
    	return false;
    }

    //check for matched pairs of parentheses
    int right = 0;
    int left = 0;
    for(int i=0; i < size; i++ )
    {
    	string test(newick,i,1);         //test is constructed with 1 character from position i
    	if(test == "(") {left++;}
    	else if (test == ")") {right++;}
    }

    //cout << "\trt = " << right << " and lf = " << left << "\n";
    if(right != left)
    {
    	cout << "\t\tdataset::validate::ERROR...parentheses are not balanced!\n";
    	return false;
    }
    return true;
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
void forestry::show_tree(int tree_id)
{
	my_forest->population[tree_id]->euler_tour();
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
string forestry::get_newick(int tree_id)
{
	return my_forest->population[tree_id]->euler_tour("newick");
}

//------------------------------------------------------------------------------
//
//------------------------------------------------------------------------------
tree* forestry::get_tree(int i)
{
	// *** add a test for i from 0 to (numtrees-1)!!!!
	// if out of range print an error!!!!
	tree* tree_ptr = my_forest->population[i];
    return tree_ptr;
}

//------------------------------------------------------------------------------
// public
//------------------------------------------------------------------------------
void forestry::check_trees()
{
	   cout << "\tCompare newick and BTM: they should be identical!\n\n";
	   for (int i=0; i < num_trees; i++)
	   {
	      cout << "\tNewick " << (i+1) << " = " << newick_trees[i] << "\n";
	      cout << "\tBTM    " << (i+1) << " = ";
	      my_forest->population[i]->euler_tour();
	      cout << "\n\n";
	   }
}

//------------------------------------------------------------------------------
// public:
//------------------------------------------------------------------------------
void forestry::delete_trees()
{
	tree* a_tree;
	//NOTE: start at ***i=1*** b/c first tree gets destroyed by destructor; see notes in tree.cpp
	for(unsigned int i=1; i<my_forest->population.size(); ++i)
	{
       a_tree = my_forest->population[i];
       my_forest->population[i] = 0;
       delete a_tree;
	}
}


