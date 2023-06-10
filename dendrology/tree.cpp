#include <iostream>
#include <stdlib.h>
#include <string>
#include <string>
#include <fstream>
#include <sstream>
#include "tree.h"
//#include "dendrology/DNA_node.h"

using namespace std;   // allows direct reference to objects such as cout

//--------------------------------------------------------
//definition of a STATIC variable in the implementation file
//--------------------------------------------------------
vector<tree*> tree::population;

//NOTE1: Population vector is defined as a static ...

//NOTE2: The first tree object is NOT abstract; a copy of itself **IS** loaded into
//       the static population when it is constructed.  This is imp when it is time to
//       delete the tree objects from memory!!!
//       -  The first tree object will get deleted when the constructor is called!
//       -  All other tree objects need to be deleted manually from the heap by a
//          "delete_trees" function b/c destroying the population vector only destroys
//          the pointers to them!!!!
//       -  The manual destruction, by using "delete" key work must therefore start at index 1
//          instead of zero for the vector, or a "double free" error is obtained

//ADD some initializer lists later (after the imp stuff is working!!!!)

//------------------------------------------------------------------------------
// CONSTRUCTOR:
//------------------------------------------------------------------------------
tree::tree():num_nodes(-1)
{
	//**********NOTES************//
	// 1: I think JRM version of code has the push_back enabled
	// 2: This might be the source of the issue that causes incompatibility between script and options interface
	//***************************//

   //population.push_back(this);
   //cout << "\n\tDEFAULT CONSTRUCTOR::tree object: size of populatin = " << population.size() << " and capacity = "<<  population.capacity()  << "\n\n";
}

//------------------------------------------------------------------------------
// CONSTRUCTOR:
//------------------------------------------------------------------------------
tree::tree(string input_str):num_nodes(-1)
{
	newick_tree = input_str;
	rooted = 0; // this will be set to 1 when the tree is read as a BTM
	// new stuff
    cout << "Loading tree: ";
    set_num_nodes();
    populate_nodes();
    cout << " ...creating BTM";
    create_btm2();
    cout << " ...traversing";
    traverse();
    cout << " ...loading branch lengths";
    load_blengths();
    cout << " ...done.\n\n";
    //set_node_heights();
    cout << "BTM = ";
    euler_tour();
    cout << "\n\n";
}

//-------------------------------------------------------
// destructor:
//-------------------------------------------------------
tree::~tree()
{
	// if(num_nodes != -1) {euler_tour();}
	//cout << "  | A tree obj. has been destroyed.\n";
}

//------------------------------------------------------------------------------
// public: some basic summary statistics for a given tree object
//------------------------------------------------------------------------------
void tree::stats()
{
	cout << "Tree Statistics:\n";
	cout << "\tnumber of OTUs = " << num_otu << "\n";
	cout << "\tnumber of nodes = " << num_nodes << "\n";
	cout << "\tnumber of branches = " << (num_otu + (num_otu-2)) << " [rooted tree only]\n";  //only for rooted tree
	cout << "\trooted = [all trees must be rooted at this time]\n";
	cout << "\troot node ID = " << get_root() << "\n";
	cout << "\ttree height = " << tree_height << "\n";    //for unrooted trees this is the height of the basal trichotomy
	cout << "\t[Add nodes listed by their height values]\n";
}


//------------------------------------------------------------------------------
// public:
//------------------------------------------------------------------------------
void tree::set_newick(string input_str)
{
	newick_tree = input_str;
}


//------------------------------------------------------------------------------
// public: [? combine with method below?]
//------------------------------------------------------------------------------
void tree::set_num_nodes()
{
	int pos = 0;
	int count = 0;
	while (newick_tree[pos] != ';')
	{
		if (newick_tree[pos] == ',') {count++;}
		pos++;
	}
	num_nodes = count + count+1;
	//cout << "\nnum nodes = " << num_nodes << "\n";
}


//------------------------------------------------------------------------------
// public:
//------------------------------------------------------------------------------
void tree::populate_nodes()
{
	// NOTE: I have abandoning derived node types!!!!

	node n;
    nodes.resize(num_nodes);
	for(int i=0; i<num_nodes; i++)
	{
		nodes.push_back(n);
		traversal_order.push_back(-1);
	}
}


//------------------------------------------------------------------------------
// public: non-recursive algorithm to construct a binary tree model (BTM)
//------------------------------------------------------------------------------
// OUTDATED !!!! [Possibly BROKEN now!!!]
// WARNING: Does not handle branch (#) or clade ($) marks !!!!
void tree::create_btm()
{
//NOTE: problems with right descendant!!!
//NOTE: this code does not work properly an is now outdated!!!!

	num_otu = 0;         // this is a global variable: initialize to 0
	the_root = -1;       // initial another global variable
	int pos=0;           // position in newick tree string
	int p=0, q=1;        // p and q index nodes in BTM
	int j=0;             // j is index of char str
	char str[100];       // temp string for taxon labels
	bool colon=false;    // flag for colon that separates OTU label from branch length

	//cout << "\ntree: " << newick_tree << "\n\n";

	while (newick_tree[pos] != ';')
	{
		if (newick_tree[pos] == '(')
		{
			//cout << "\n\"(\"\n";
			//cout << "p: " << p << "\n";
			//cout << "q: " << q << "\n";

			nodes[p].set_label("internal");
			//cout << "The label of node " << p << " is now INTERNAL.\n";

			nodes[q].set_anc(p);     // set ancestor of q to p
			nodes[p].set_left(q);    // set left descendant of p to q

			//int x = nodes[q].get_anc();
			//cout << "anc of node "  << q << " is " << x << "\n";
			//int y = nodes[p].get_left();
			//cout << "left of node "  << p << " is " << y << "\n";

			nodes[p].set_id(p);
			p++;
			q++;
			pos++;

		}
		else if (newick_tree[pos] == ',')
		{
			//cout << "\n\",\"\n";
			//cout << "p: " << p << "\n";
			//cout << "q: " << q << "\n";

			/***    OTU label bookkeeping    ***/
			str[j] = '\0';                      // add a ???? to the string
			j=0;                                // done with OTU label, so reset its index
			nodes[p].set_label(str);            // set label of node
			//cout << "label of " << p << " is " << str << "\n";
			//cout << "label of node " << p << " is now " << nodes[p].get_label() << "\n\n";

			/***    misc bookkeeping    ***/
			num_otu++;                     // increment number of OTUs
			colon = false;                 // set/reset flag for colon to false so OTU labels will get read properly below

			int k = nodes[p].get_anc();    // k is ancestor of p
			int l = nodes[k].get_left();   // l is left descendant of k

			/// PART A ///
			if(p==l)                       // if p is the left descendant of its ancestor
			{
				//cout << "p == l: " << p << " == " << l << "\n";
				//cout << "k = " << k << "\n";

				int k = nodes[p].get_anc(); // k is ancestor of current p
				nodes[k].set_right(q);      // set right descendant of ancestor of p to q: GOOD
				nodes[q].set_anc(k);        // set ancestor of q to ancestor of p: GOOD

				//int x = nodes[k].get_right();
				//cout << "right of node "  << k << " is " << x << "\n";
				//int y = nodes[q].get_anc();
				//cout << "anc of node "  << q << " is " << y << "\n";
			}
			/// PART B ///
			else                            // if p is right descendant of its ancestor
			{
				// BUG MUST BE IN HERE?
				//cout << "p != l: " << p << " != " << l << "\n";
				int k = nodes[p].get_anc();  // k is ancestor of p
				int m = nodes[k].get_anc();  // m is ancestor of ancestor of p (i.e., k)
				//nodes[m].set_right(q);       //

				nodes[k].set_right(p);
				nodes[p].set_anc(k);
				// this part causes problem with central root in tree: put within an if block to fix?

				if( nodes[m].get_right() != k )
				{
					nodes[m].set_right(q);       // set right descendant of ancestor of ancestor of p to q: bad?
					nodes[q].set_anc(m);         // set ancestor of p (q) to ancestor of ancestor of p (m): good?
				}

			}
			//nodes[q].set_anc(p);
			nodes[p].set_id(p);
			p++;
			q++;
			pos++;
		}
		else if (newick_tree[pos] == ')')
		{
			//cout << "\n\")\"\n";
			//cout << "p: " << p << "\n";
			//cout << "q: " << q << "\n";
			// TRY: set p to ancestor of q
			colon = false;   // reset flag to false
			pos++;
		}
		else if (newick_tree[pos] == ':')
		{
			//cout << "\n\":\"\n";
			//cout << "p: " << p << "\n";
			//cout << "q: " << q << "\n";
			colon = true;  // set flag to true when a colon is observed
			pos++;
		}
		else
		{

			if(colon == false) // if a colon has not been seen yet, keep adding characters to the taxon label
			{
				str[j] = newick_tree[pos]; j++;  // reading label of p
			}
			pos++;
		}
		/// end if statements
	}
	/// end while loop

	str[j] = '\0';
	nodes[p].set_label(str);
	nodes[p].set_id(p);
	num_otu++;
	//cout << "label: " << str << "\n\n";
	//cout << "The number of OTUs in the tree is " << num_otu << "\n";
	num_nodes = q;
	rooted = 1;
}


//------------------------------------------------------------------------------
// PUBLIC: alternative non-recursive algorithm to construct a binary tree model (BTM)
// NOTE: This ONLY works on a rooted (fully bifurcating) newick tree!!!!
//------------------------------------------------------------------------------
void tree::create_btm2()
{
	num_otu = 0;         // this is a global variable: initialize to 0
	the_root = -1;       // initial another global variable
	int pos=0;           // position in newick tree string
	int p=0, q=1;        // p and q index nodes in BTM
	int j=0;             // j is index of char str
	char str[100];       // temp string for taxon labels
	bool colon=false;    // flag for colon that separates OTU label from branch length

	while (newick_tree[pos] != ';')
	{
		if (newick_tree[pos] == '(')
		{
			nodes[p].set_label("internal");
			nodes[p].set_id(p);
			nodes[q].set_anc(p);     // set ancestor of q to p
			nodes[p].set_left(q);    // set left descendant of p to q

			int right = get_right_descendant(p, newick_tree);
            // JPB::17Nov2015: If at the ti of the tree, the above method will
            // return a -1.  In this case there is no right descendant for node p,
            // and thus there will be no node[right] that needs an ancestor set.
            // Solution was to onlu update below if right != -1.
            if(right != -1)
            {
                nodes[p].set_right(right);
                nodes[right].set_anc(p);
            }
			p++;
			q++;
			pos++;
		}
		else if (newick_tree[pos] == ',')  // working on tip nodes in here
		{
			/***    OTU label bookkeeping    ***/
			str[j] = '\0';                      // add a ... to the string
			j=0;                                // done with OTU label, so reset its index
			nodes[p].set_label(str);            // set label of node

			/***    misc bookkeeping    ***/
			num_otu++;                     // increment number of OTUs
			colon = false;                 // set/reset flag for colon to false so OTU labels will get read properly below

			int k = nodes[p].get_anc();    // k is ancestor of p
			int l = nodes[k].get_left();   // l is left descendant of k

			/// PART A: when p is the left descendant of its ancestor
			if(p==l)
			{
				int k = nodes[p].get_anc(); // k is ancestor of current p
				nodes[k].set_right(q);      // set right descendant of ancestor of p to q: GOOD
				nodes[q].set_anc(k);        // set ancestor of q to ancestor of p: GOOD
			}
			/// PART B: when p is right descendant of its ancestor
			else
			{
				int k = nodes[p].get_anc();  // k is ancestor of p
				int m = nodes[k].get_anc();  // m is ancestor of ancestor of p (i.e., k)
				nodes[k].set_right(p);
				nodes[p].set_anc(k);
			}
			nodes[p].set_id(p);
			p++;
			q++;
			pos++;
		}
		else if (newick_tree[pos] == ')')
		{
			colon = false;   // reset flag to false
			pos++;
		}
		else if (newick_tree[pos] == ':' || newick_tree[pos] == '#'   ||  newick_tree[pos] == '$' )
		{
			colon = true;  // set flag to true when a colon, or some other mark, is observed
			pos++;
		}
		else
		{
			if(colon == false) // if a colon has not been seen yet, keep adding characters to the taxon label
			{
				str[j] = newick_tree[pos]; j++;  // reading label of p
			}
			pos++;
		}
		/// end if statements

	}
	/// end while loop

	str[j] = '\0';
	nodes[p].set_label(str);
	nodes[p].set_id(p);
	num_otu++;
	//cout << "label: " << str << "\n\n";
	//cout << "The number of OTUs in the tree is " << num_otu << "\n";
	num_nodes = q;
	rooted = 1;
}


//------------------------------------------------------------------------------
// public:
//------------------------------------------------------------------------------
int tree::get_right_descendant(int q, string local_tree)
{
	//NOTE: creates array each time; if works, then fix so that this is done only once.
	int pos=0;
	int p=0;
	bool tip=false;


	int size = local_tree.size();
	int node_index[size];
	int rt_descend[size];

	while (local_tree[pos] != '\0')
	{
		if (local_tree[pos] == '(')
		{
			node_index[pos] = p;
			p++;
			pos++;
		}
		else if (local_tree[pos] == ',')
		{
			node_index[pos] = p;
			p++;
			pos++;
		}
		else if (local_tree[pos] == ';')
		{
			node_index[pos] = p;
			p++;
			pos++;
		}
		else
		{
			node_index[pos] = -1;
			pos++;
		}
	}

	pos = 0;

	/*
	  cout << "\nSize of newick string is " << size << "\n";
	  cout << "\nnode indices.....\n";
	  for(int i=0; i<size; i++)
	  {
	  cout << "At pos = " << i << " the node ID is " << node_index[i] << "\n";
	  }
	  cout << "\n";
	*/

	// NOTE: this part could be modified so that tip node differs from rt descend==tip
	//
	// First find position of node q within the newick tree sting
	// Second add 1 (pos+1) because that will be start of descendant clade
	for(int i=0; i<size; i++)
	{
		if (node_index[i] == q)
		{
			pos = i;
			pos = pos+1;
			//cout << "\ngiven " << q << " I am looking at pos " << pos << "\n";
			if (local_tree[pos] != '(')
			{
				//cout << "\npos = " << pos << "\n";
				//cout << "char = " << newick_tree[pos] << "\n";
				//cout << "\n\nRight descendant of node " << q << " is a tip.\n";
				tip = true;
			}
			break;
		}
	}

	// Check if descendant clade is a tip (i.e., starts with "(")
	// If it is NOT a tip, then there is no need to do anything more
	// If next descendant it not a "(" then the clade is actually
	// a tip and the function will return a "-1".
	if (tip == true) {return -1;}

	// Look for the position of the rt descendant within the descendant clade
	// This starts on the next left bracket of the descendant clade
	int counter = 0;
	int look_here = 0;
	bool flag = false;
	for(int i=pos; i<size; pos++)
	{
		// start by summing rt an lt brackets as we move to the rt
		if (flag == false)
		{
			if (local_tree[pos] == '(') {counter ++;}
			else if (local_tree[pos] == ')') {counter --;}
			else {/*no nothing*/}
		}

		// when rt = lf the counter == 0
		// set flag true
		if (counter == 0)
		{
			flag = true;
			//cout << "Found something interesting\n";
		}

		// now find position of first "," after the left most ")"
		// of the descendant clade
		//if (flag == true && (newick_tree[pos] == '(' || newick_tree[pos] == ';'))
		//****chnage here: check for bugs******
		if (flag == true && (local_tree[pos] == ',' || local_tree[pos] == ';'))
		{
			look_here = pos;  // look at pos of fist ',' after counter==0
			break;            // and stop the search
		}
	}

	//cout << "Right descendant of " << q << " is node " << node_index[look_here] << "\n";
	return node_index[look_here]+1;
}


//------------------------------------------------------------------------------
// public: load branch lengths to nodes in a fully specified BTM
//------------------------------------------------------------------------------
void tree::load_blengths()
{
	int pos=0;           // position in newick tree string
	int i=0;             // index node_ID numbers in the traversal_order (DOWNPASS) array
	int k=0;             // h is index for char str2
	char str[30];        // temp string for branch lengths
	bool colon=false;    // flag for colon that separates OTU label from branch length

	//cout << "\n";
	while (newick_tree[pos] != ';')
	{
		//cout << newick_tree[pos];
		//cout << pos << " ";
		if (newick_tree[pos] == '(')
		{
			pos++;
		}
		else if (newick_tree[pos] == ',')
		{
			/***  branch length bookkeeping  ***/
			str[k] = '\0';                            // put \0 on end of current str
			//cout << "\n ,blength =" << str << "\n";
			k=0;                                      // reset k back to 0
			double temp_double = -1;                  // default value = -1
			istringstream my_iss(str);                // open a string stream for the float copied to str
			my_iss >> temp_double;                    // assign str to an *float* type variable
			int node_id = get_downpass_node(i);       // NEW: to safely work with the traversal order vector
			//cout << "Node(" << node_id << "): Branch length = " << temp_double << "\n";
			nodes[node_id].set_blength(temp_double);  // set function used to set float as a node branch length
			colon = false;                            // set/reset flag for colon to false so OTU labels will get read properly below
			pos++;
			i++;
		}
		else if (newick_tree[pos] == ')')
		{
			colon = false;                            //
			str[k] = '\0';                            //
			k=0;                                      //
			double temp_double = -1;                  // default value = -1
			istringstream my_iss(str);                // open a string stream for the float copied to str2
			my_iss >> temp_double;                    // assign str2 to an *float* type variable
			int node_id = get_downpass_node(i);       // NEW: to safely work with the traversal order vector
			//cout << "Node(" << node_id << "): Branch length = " << temp_double << "\n";
			nodes[node_id].set_blength(temp_double);  // set function used to set float as a node branch length
			//cout << "\n )blength =" << str << "\n";
			pos++;
			i++;
		}
		else if (newick_tree[pos] == ':')
		{
			colon = true;  // set flag to true when a colon is observed
			pos++;
		}
		else
		{
			if (colon == true)  // collect everything after the colon in a str
			{
				str[k] = newick_tree[pos];
				k++;
			}
			pos++;
		}
		/// end if statements
	}
	/// end while loop
}

//------------------------------------------------------------------------------
// public: load branch marks identified within a newick tree by the "#" symbol
//------------------------------------------------------------------------------
int tree::load_branch_marks()
{
	//NOTE1: Any number of branch marks will be allowed in a given tree
	//NOTE2: When branch marks are read from a different tree than the primary one (having
	//       the branch lengths, then the two newick trees must be check for identity! Such
	//       a check is NOT done here!!!

	int pos=0;           // position in newick tree string
	int i=0;             // index node_ID numbers in the traversal_order (DOWNPASS) array
	int k=0;             // h is index for char str2
	int counter = 0;     // count number of instances of the "#" symbol in the tree file
	char str[30];        // temp string for branch lengths
	bool hash=false;     // flag for hash-symbol that separates OTU label from branch mark

	//cout << "\n";
	while (newick_tree[pos] != ';')
	{
		//cout << newick_tree[pos];
		//cout << pos << " ";
		if (newick_tree[pos] == '(')
		{
			pos++;
		}
		else if (newick_tree[pos] == ',')
		{
			/***  branch mark bookkeeping  ***/
			str[k] = '\0';                               // put \0 on end of current str
			//cout << "\n ,branch-mark =" << str << "\n";
			k=0;                                         // reset k back to 0
			int temp_int = 0;                            // default value = 0
			istringstream my_iss(str);                   // open a string stream for the float copied to str
			my_iss >> temp_int;                          // assign str to an *int* type variable
			int node_id = get_downpass_node(i);          // to safely work with the traversal order vector
			//cout << "Node(" << node_id << "): Branch mark = " << temp_int << "\n";
			nodes[node_id].set_branch_class(temp_int);   // set function used to set float as a node branch length
			hash = false;                                // set/reset flag for hash to false so OTU labels will get read properly below
			pos++;
			i++;
		}
		else if (newick_tree[pos] == ')')
		{
			hash = false;                               //
			str[k] = '\0';                              //
			k=0;                                        //
			int temp_int = 0;                           // default value = 0
			istringstream my_iss(str);                  // open a string stream for the float copied to str2
			my_iss >> temp_int;                         // assign str2 to an *float* type variable
			int node_id = get_downpass_node(i);         // to safely work with the traversal order vector
			//cout << "Node(" << node_id << "): Branch mark = " << temp_int << "\n";
			nodes[node_id].set_branch_class(temp_int);  // set function used to set float as a node branch length
			//cout << "\n )branch-mark =" << str << "\n";
			pos++;
			i++;
		}
		else if (newick_tree[pos] == '#')
		{
			hash = true;  // set flag to true when a colon is observed
			counter++;
			pos++;
		}
		else
		{
			if (hash == true)
			{
				str[k] = newick_tree[pos];
				k++;
			}
			pos++;
		}
		/// end if statements
	}
	/// end while loop

	if(counter == 0)      return 0;
	else if (counter > 0) return 1;
	else                  return -1;
}


//------------------------------------------------------------------------------
// public: load clade marks identified within a newick tree by the "$" symbol
//------------------------------------------------------------------------------
int tree::load_clade_marks()
{
	//NOTE1: Currently only ONE clade marks will be allowed in a given tree
	//NOTE2: When clade marks are read from a different tree than the primary one (having
	//       the branch lengths, then the two newick trees must be checked for identity! Such
	//       a check is NOT done here!!!

	int pos=0;           // position in newick tree string
	int i=0;             // index node_ID numbers in the traversal_order (DOWNPASS) array
	int k=0;             // h is index for char str2
	int counter = 0;     // counter for number of clade marks ($ symbols in the tree); can be zero or 1 only
	char str[30];        // temp string for branch lengths
	bool dollar=false;   // flag for hash-symbol that separates OTU label from branch mark

	//cout << "\n";
	while (newick_tree[pos] != ';')
	{
		//cout << newick_tree[pos];
		//cout << pos << " ";
		if (newick_tree[pos] == '(')
		{
			pos++;
		}
		else if (newick_tree[pos] == ',')
		{
			/***  clade-mark bookkeeping  ***/
			str[k] = '\0';                               // put \0 on end of current str
			//cout << "\n , clade-mark =" << str << "\n";
			k=0;                                         // reset k back to 0
			int temp_int = 0;                            // default value = 0
			istringstream my_iss(str);                   // open a string stream for the float copied to str
			my_iss >> temp_int;                          // assign str to an *int* type variable
			int node_id = get_downpass_node(i);          // to safely work with the traversal order vector
			//cout << "Node(" << node_id << "): Clade mark = " << temp_int << "\n";
			nodes[node_id].set_clade_class(temp_int);    // set function used to set float as a node branch length
			dollar = false;                              // set/reset flag for hash to false so OTU labels will get read properly below
			pos++;
			i++;
		}
		else if (newick_tree[pos] == ')')
		{
			dollar = false;                             //
			str[k] = '\0';                              //
			k=0;                                        //
			int temp_int = 0;                           // default value = 0
			istringstream my_iss(str);                  // open a string stream for the float copied to str2
			my_iss >> temp_int;                         // assign str2 to an *float* type variable
			int node_id = get_downpass_node(i);         // to safely work with the traversal order vector
			//cout << "Node(" << node_id << "): Clade mark = " << temp_int << "\n";
			nodes[node_id].set_clade_class(temp_int);   // set function used to set float as a node branch length
			//cout << "\n )clade =" << str << "\n";
			pos++;
			i++;
		}
		else if (newick_tree[pos] == '$')
		{
			dollar = true;  // set flag to true when a colon is observed
			counter ++;
			pos++;
		}
		else
		{
			if (dollar == true)
			{
				str[k] = newick_tree[pos];
				k++;
			}
			pos++;
		}
		/// end if statements
	}
	/// end while loop


	//check that there is only on dollar sign in the tree
	if(counter > 1)
	{
		cout << "ERROR::tree::load_clade_marks: Only 1 clade mark permitted at this time!\n";
		cout << "Sorry, this program must exit.\n";
		exit(0);
	}
	else if (counter == 0)
	{
		cout << "CAUTION::tree::load_clade_marks: Expecting a clade mark and did not find one!\n";
		cout << "Attempting to proceed anyway(you have been warned)!\n";
	}


	// save node ID having the (assumed) single dollar-symbol
	int clade_class = 0;
	int foreground_node = -1;
	for(int i=(num_nodes-1); i>=0; i--)
	{
		clade_class = nodes[i].get_clade_class();
		if(clade_class != 0) foreground_node = i;
	}

	//change clade label according to the single dollar-symbol
	int ancestor_class = 0;
	int right = -1;
	int left = -1;
	int node;
	for(int i=(num_nodes-1); i>=0; i--)
	{
		node = get_downpass_node(i);
		//cout << "working on node = " << node << "\n";
		ancestor_class = nodes[node].get_clade_class();
		//cout << "class label = " << ancestor_class << "\n";
		right = nodes[node].get_right();
		left = nodes[node].get_left();
		//*** Only allow one mark; so, regardless of the int, set the mark = 1 ***//
		if(ancestor_class != 0)
		{
			nodes[node].set_clade_class(1);
			if (right != -1) nodes[right].set_clade_class(1);
			if (left != -1)  nodes[left].set_clade_class(1);
		}
	}

	return foreground_node;
}


//------------------------------------------------------------------------------
// public: shows all the nodes; be careful b/c this is way too much for a big tree!!!
//------------------------------------------------------------------------------
void tree::show_nodes()
{
	cout << "\nNewick      : " << newick_tree << "\n";
	cout << "\nThe BTM is  : ";
	euler_tour();
	cout << "\n\n";
	cout << "\nThe number of OTUs is: " << num_otu << "\n";
	cout << "\nThe total number of nodes (including tips) is: " << num_nodes << "\n";

	for(int i=0; i<num_nodes; i++)
	{
		cout << "\n\n NODE NUMBER " << i << "\n";
		nodes[i].show();
		//if (nodes[i].get_anc() == -1) {the_root = i;}
	}
}


//------------------------------------------------------------------------------
// public: Traverse the tree by calling "uppass"; note that "pass" is recursive
//------------------------------------------------------------------------------
void tree::traverse()
{
	index=0;
	//cout << "Starting the traversal at node 0: \n";
	uppass(nodes[0]);  // nodes[0] should always be the root
	// start traversal of tree at root
}

//------------------------------------------------------------------------------
// private: recursive function; if a node is NOT a tip, it calls itself for LF & RT descendants
//------------------------------------------------------------------------------
void tree::uppass(node &input_node)
{
	bool tip;
	int id = input_node.get_id();

	//cout << "\nuppass: \n";
	if( (input_node.get_right() == -1) && (input_node.get_left() == -1) )
	{
		//cout << "Node " << id << " is a tip node\n";
		//cout << "Label of this tip is " << input_node.get_label() << "\n";
	}
	else
	{
		//cout << "Node " << id << " is NOT a tip node\n";
		int left = input_node.get_left();
		int right = input_node.get_right();
		//cout << "The left descendant is node " << left << "\n";
		//cout << "The right descendant is node " << right << "\n";
		uppass( nodes[left] );
		uppass( nodes[right] );
	}
	//cout << "index: " << index << "\n";
	traversal_order[index] = id;
	index++;
}

//------------------------------------------------------------------------------
// public: show the DOWNPASS traversal order as given in "traversal_order" array
//------------------------------------------------------------------------------
void tree::print_downpass()
{
	//cout << "\nSize of traversal vector = " << traversal_order.size() << "\n";
	cout << "\ntraversal order: ";
	for(int i=0; i<num_nodes; i++) { cout << traversal_order[i] << ", "; }
	cout << "\n";


	cout << "\nORDER OF NODES IN A DOWNPASS\n\n";
	for(int i=0; i<num_nodes; i++)
	{
		cout << i << " is node " << traversal_order[i] << ", label = " << nodes[traversal_order[i]].get_label();
		cout << ", \theight = " << nodes[traversal_order[i]].get_height() << "\n";
	}
}

//------------------------------------------------------------------------------
// public: show the UPPASS traversal order as given in "traversal_order" array
//------------------------------------------------------------------------------
void tree::print_uppass()
{
	cout << "\nORDER OF NODES IN A UP-PASS\n\n";
	for(int i=(num_nodes-1); i>=0; i--)
	{
		cout << i << " is node " << traversal_order[i] << ", label = " << nodes[traversal_order[i]].get_label() << "\n";
	}
}


//------------------------------------------------------------------------------
// public:
//------------------------------------------------------------------------------
void tree::show_heights()
{
	cout << "\tnode heights:  ";
	for(int i=0; i<num_nodes; i++)
	{
		cout << nodes[traversal_order[i]].get_height() << "(N" << traversal_order[i] << "), ";
	}
	cout << "\n";
}


//------------------------------------------------------------------------------
// public: show to screen a variety of info related to branch classes
//------------------------------------------------------------------------------
void tree::show_branch_classes()
{
	// show a newick tree with IDs for branch classes
	cout << "\n\nBranch class details:\n";
	cout << "\tNewick tree with branches labeled with class IDs of a \"branch class model\":\n";
	cout << "\t[node]:\t\t";
	labeled_euler_tour(1);
	cout << "\n\t<branch class>:\t";
	labeled_euler_tour(0);
	cout << "\n\t{clade class}:\t";
	labeled_euler_tour(2);
	cout << "\n\n";

	// show a node list with info relevant to setting the branch class IDs
	bool tip = false;
	cout << "\tNode  \tLeft  \tRight  \tBranch-Class \tClade-class \n";
	for(int i=0; i<num_nodes; i++)
	{
		cout << "\t" << nodes[i].get_id();  if (nodes[i].get_anc() == -1) {cout << "*";}
		if(nodes[i].is_tip()) {tip = true;}
		if(tip)
		{
			cout << "\t(Tip)";
			cout << "\t(Tip)";
		}
		else
		{
			cout << "\t" << nodes[i].get_left();
			cout << "\t" << nodes[i].get_right();
		}
		cout << "\t" << nodes[i].get_branch_class();
		cout << "\t\t" << nodes[i].get_clade_class();
		if (nodes[i].get_anc() == -1) {cout << "\t\t(*This is the root node)";}
		cout << "\n";
	}
	cout << "\n";

}

//------------------------------------------------------------------------------
// public:
//------------------------------------------------------------------------------
string tree::get_node_label(int node_num)
{
	string label = nodes[node_num].get_label();
	return label;
}

//------------------------------------------------------------------------------
// public:
//------------------------------------------------------------------------------
int tree::get_num_nodes()
{
	return num_nodes;
}

//------------------------------------------------------------------------------
// public:
//------------------------------------------------------------------------------
int tree::get_num_otu()
{
	return num_otu;
}

//------------------------------------------------------------------------------
// public:
//------------------------------------------------------------------------------
int tree::get_root()
{
	if(the_root == -1)
	{
		for(int i=0; i<num_nodes; i++)
		{
			if (nodes[i].get_anc() == -1) {the_root = i;}
		}
	}
	//cout<< "\nThe root node is " << the_root << "\n\n";
	return the_root;
}


//------------------------------------------------------------------------------
// public:
//------------------------------------------------------------------------------
int tree::get_downpass_node(int rank)
{
	return traversal_order[rank];
}


//------------------------------------------------------------------------------
// public: Euler tour algorithm
//------------------------------------------------------------------------------
void tree::euler_tour()
{
	// Euler tour starts at root.
	// cout << "\nStarting Euler tour....\n";
	int i = get_root();
	// cout << "Starting algorithm at node " << i << "\n\n";
	euler(i);
	// cout << ";\n\n" << "Finished Euler tree tour!!!!\n";
}


//------------------------------------------------------------------------------
// private: the is the RECURSIVE part of the Euler tour algorithm
//------------------------------------------------------------------------------
// NOTE: when tree has no blengths the default value is set to -1; hence,
// trees that have blength != -1 have a branch length of interest.
void tree::euler(int i)
{
	if(nodes[i].is_tip())
	{
		cout << nodes[i].get_label();
		// if tip node has branch lengths, then print them
		if( nodes[i].get_blength() != -1 )
		{
			cout << ":" << nodes[i].get_blength();
		}
	}
	else
	{
		cout << "(";                  // "on the left of node i"
		euler(nodes[i].get_left());   // recursion to the left
		cout << ",";                  // "from below node i"
		euler(nodes[i].get_right());  // recursion to the right
		cout << ")";                  // "on the right of node i"

		// if internal node is not the root AND it has branch lengths, then print
		if(i != get_root() && nodes[i].get_blength() != -1)
		{
			cout << ":" << nodes[i].get_blength();
		}
	}
}


//------------------------------------------------------------------------------
// public: Euler tour algorithm for newick tree with BRANCH CLASS IDs
//------------------------------------------------------------------------------
void tree::labeled_euler_tour(int flag)
{
	int i = get_root();
	labeled_euler(i, flag);
}


//------------------------------------------------------------------------------
// private: The is the RECURSIVE part of the labeled_Euler tour algorithm
//------------------------------------------------------------------------------
void tree::labeled_euler(int i, int flag)
{
	if(nodes[i].is_tip())
	{
		cout << nodes[i].get_label();
		if(flag == 0)      { cout << ":<" << nodes[i].get_branch_class() << ">"; }
		else if(flag ==1)  { cout << ":[" << nodes[i].get_id() << "]"; }
		else if(flag ==2)  { cout << ":{" << nodes[i].get_clade_class() << "}"; }
		else               { cout << ":" << nodes[i].get_id() << "/" << nodes[i].get_branch_class(); }
	}
	else
	{
		cout << "(";                  // "on the left of node i"
		labeled_euler(nodes[i].get_left(), flag);   // recursion to the left
		cout << ",";                  // "from below node i"
		labeled_euler(nodes[i].get_right(), flag);  // recursion to the right
		cout << ")";                  // "on the right of node i"
		if(i != get_root() )
		{
			if(flag == 0)      { cout << ":<" << nodes[i].get_branch_class() << ">"; }
			else if(flag ==1)  { cout << ":[" << nodes[i].get_id() << "]"; }
			else if(flag ==2)  { cout << ":{" << nodes[i].get_clade_class() << "}"; }
			else               { cout << ":" << nodes[i].get_id() << "/" << nodes[i].get_branch_class(); }
		}
	}
}


//------------------------------------------------------------------------------
// public: Overloaded Euler tour
//------------------------------------------------------------------------------
string tree::euler_tour(string str)
{
	string str2;
	string newick_str;
	// Euler tour starts at root.
	int i = get_root();
	newick_str = euler2(i, str2);
	return newick_str;
}

//------------------------------------------------------------------------------
// PIRIVATE: the is the RECURSIVE part of the Euler tour algorithm
//------------------------------------------------------------------------------
// NOTE: This function differs in that it generates a string that contains
//       trees in newick format.
// NOTE: this is called by overloaded euler_tour() function

string tree::euler2(int i, string str)
{
	// need to open and use a string stream to handle the blenghts
	// b/c blengths are floats, and they must be put into a string
	stringstream ss;

	if(nodes[i].is_tip())
	{
		str += nodes[i].get_label();
		// if tip node has branch lengths, then print them
		if( nodes[i].get_blength() != -1 )
		{
			str += ":";
			ss << nodes[i].get_blength();  // use string stream to handle the float here
			str += ss.str();               // access it as a string using the .str() function
			ss.str("");                    // empty the string
		}
	}
	else
	{
		str += "(";                               // "on the left of node i"
		str = euler2(nodes[i].get_left(), str);   // recursion to the left
		str += ",";                               // "from below node i"
		str = euler2(nodes[i].get_right(), str);  // recursion to the right
		str +=  ")";                              // "on the right of node i"

		// if internal node is not the root AND it has branch lengths, then print
		if(i != get_root() && nodes[i].get_blength() != -1)
		{
			str += ":";
			ss << nodes[i].get_blength();  // use string stream to handle the float here
			str += ss.str();               // access it as a string using the .str() function
			ss.str("");                    // empty the string
		}
	}
	return str;
}



//------------------------------------------------------------------------------
// PRIVATE:
//------------------------------------------------------------------------------
void tree::set_node_heights()
{
	//cout << "\tnode heights:  ";
	tree_height = -1;
	for(int i=0; i<num_nodes; i++)
	{
		if(nodes[i].is_tip())
		{
			nodes[i].set_height(0);
			recursive_heights(0, nodes[i]);
		}
	}
}

//------------------------------------------------------------------------------
// PRIVATE: recursive function to assign height values to all nodes of the tree
//------------------------------------------------------------------------------
void tree::recursive_heights(int ht, node &input_node)
{
	int current_ht = input_node.get_height();

	if(ht > current_ht)
	{
		input_node.set_height(ht);
		//cout << " " << ht << " ";
		//cout << " " << ht << "(" << input_node.get_id() << ")  ";
		if(ht > tree_height)
		{
			tree_height = ht;
		}
	}

	if(input_node.is_root()) {/*do nothing*/}
	else
	{
		ht++;
		int anc = input_node.get_anc();
		recursive_heights(ht, nodes[anc]);
	}
}


//------------------------------------------------------------------------------
// public:
//------------------------------------------------------------------------------
void tree::delete_trees()
{
	tree* a_tree;
	//Note: Start at i=1 b/c first tree gets destroyed by destructor; see notes in tree.cpp
    //      Thus, nothing should happen here when there is only one tree in the tree population vector.
	for(unsigned int i=1; i<population.size(); ++i)
	{
		a_tree = population[i];   // this line is generating the problem!!!
		population[i] = 0;
		delete a_tree;
	}
}

//------------------------------------------------------------------------------
// public:
//------------------------------------------------------------------------------
void tree::show_states(bool internal)
{
	cout << "Character states at node:\n";
	for(int i=0; i<num_nodes; i++)
	{
		if(nodes[i].is_tip())
		{
			cout << nodes[i].get_label() << ": " << nodes[i].get_state() << "\n";
		}
		else if(internal)
		{
			cout << "Internal node (ID " << nodes[i].get_id() << "): " << nodes[i].get_state() << "\n";
		}
	}
}

//------------------------------------------------------------------------------
// PUBLIC: Tree stays as BTM (NOT set to basal trichotomy); 1 root blength set=0
//------------------------------------------------------------------------------
void tree::deroot()
{
	// NOTES: 1. Root node actually stays the same!
	//        2. blenght to left_of_root always set = 0
	//        3. Bool variable rooted set to false

	// Check if tree is already unrooted
	if (!rooted)
	{
		cout << "CAUTION::tree::deroot: Cannot \"de-root\" an unrooted tree: no further action taken.\n";
		return;
	}

	// set up req'd variables
	double new_blength = -1;
	int left_of_root = nodes[the_root].get_left();
	int right_of_root = nodes[the_root].get_right();

	// If the tree has branch lengths (anything other than "-1") then
	// set the new branch between "left_of_root" an "right_of_root" equal
	// to the sum of the left and right branch lengths to the "old root"
	if(nodes[left_of_root].get_blength() != -1)
	{
		new_blength = nodes[left_of_root].get_blength() + nodes[right_of_root].get_blength();
	}

	nodes[right_of_root].set_blength(new_blength);  // branch to right of root gets the sum (or -1)
	nodes[left_of_root].set_blength(0.0);           // branch to left of root always gets FIXED to zero
	rooted = false;                                 // False signifies that branch to left of root is FIXED = 0

	cout << "\tTree has been de-rooted.\n";
	cout << "\t"; euler_tour(); cout << "\n";
}

//------------------------------------------------------------------------------
// PUBLIC:
//------------------------------------------------------------------------------
bool tree::is_rooted()
{
	if(rooted) {return true;}
	else       {return false;}
}

//------------------------------------------------------------------------------
// PUBLIC: Seems to work; blengths seem to be dealt with!!!
//------------------------------------------------------------------------------
void tree::set_root(int target, bool blengths, bool bmarks)
{
	//NOTES: 1. complex and need more testing
	//       2. roots along one branch; future can take a set of taxa and find the
	//          branch they represent and then call this function to root the tree!!!
	//
	// NOTE: If bmarks flag is set = "true" then any bmarks in the tree will be adjusted
	//       to suit the new location of the root.  BUT, if bmarks = "true" then this
	//       method first checks that any marks to right and let of root are equal. If
	//       such marks are not equal then the program is terminated. I.e., you cannot
	//       reroot a tree with rt and left descendants having different evolution!
	//
	//
	// NOTE perhaps code a new method that takes a pair of nodes as arguments and roots tree between them!!!


	//cout << "attempting to root tree at node " << target << ".\n";

	//bool bmarks=false;

	//----------------------------------------------------------------------------//
	//----                check for a valid target node                       ----//
	//----------------------------------------------------------------------------//
	if( target < 0  ||  target > (num_nodes-1) )
	{
		cout << "\nERROR::tree::set_root: Invalid target node(" << target << ").\n";
		cout << "- A valid target node ID must be within 0 and " << (num_nodes-1) << ".\n";
		cout << "- Proceeding with tree as currently rooted.\n";
		return;
	}

	//----------------------------------------------------------------------------//
	//---- check that the tree is not already rooted at under the target node ----//
	//----------------------------------------------------------------------------//
	int anc = nodes[target].get_anc();
	if(anc == -1)
	{
		cout << "\nTree is already rooted relative the supplied target node: no further action taken.\n";
		return;
	}

	//----------------------------------------------------------------------------//
	//----     check that user is not trying to root along a root branch      ----//
	//----------------------------------------------------------------------------//
	if(anc == 0)
	{
		cout << "\nCAUTION::tree:set_root: Attempted to place root along branch\n "
                "connected to the root node. This has no effect.\n"; /*Rotating\n"
																	   "left and right descendants of root node: Is this what you wanted?\n";*/

		//NOTE: Turned off (for now) so that fast bl optimization will work properly.
		//int new_left = nodes[the_root].get_right();
		//nodes[the_root].set_right(nodes[the_root].get_left());
		//nodes[the_root].set_left(new_left);
		return;
	}


	//----------------------------------------------------------------------------//
	//--- check that user is not trying to re-root a tree with different marks ---//
	//----------------------------------------------------------------------------//
	if(bmarks)
	{
		int rt, lf;
		//check that marks are same on both sides of the root node
		rt = nodes[0].get_right();  //cout << "rt = " << rt << "\n";
		lf = nodes[0].get_left();   //cout << "lf = " << lf << "\n";
		if(nodes[rt].get_branch_class() != nodes[lf].get_branch_class())
		{
			cout << "ERROR::tree::set_root: re-setting root in not allowed when branch mark to right of root is not equal to branch mark to left of root!\n";
			cout << "Sorry, this progam must terminate.\n";
			exit(0);
		}
	}

	//------------------------------------------------------------------//
	//---- Find the chain of nodes between target and old root node ----//
	//------------------------------------------------------------------//
	// NOTE: The "node-chain" is the chain of nodes that will/have been
	//       affected by changing the root; all others will be "OK" (I hope).
	//       Nodes affected are those that lie along direct path in tree
	//       from new root location to the old root location. Pointers affected
	//       within those nodes are the ones that point to node in this path (their
	//       direction (anc vs. descendant) must be reversed.

	vector<int> node_chain;
	int id = nodes[target].get_id();
	node_chain.push_back(id);
	while (nodes[id].get_id() != the_root) // this gets "node-chain" by working back from target node to old root
	{
		id = nodes[id].get_anc();
		node_chain.push_back(id);
	}
	node_chain.pop_back();                              // take the root node off the end of the chain, as it will be deleted from its current location in the BTM
	node_chain.insert(node_chain.begin() +1, the_root); //put the root in as the ancestor of the target node (i.e., this is the "new" location for the root!!!

	// I think the last one in the chain is OK, so the swap algorithm below will not make change on the last one [I hope]
	// This puts the "last one" (node) on the chain
	int right = nodes[the_root].get_right();
	int left = nodes[the_root].get_left();
	int size = node_chain.size();
	if(node_chain[size-1] == right) {node_chain.push_back(left);}
	else                            {node_chain.push_back(right);}

	/*
	  for(unsigned int i=0; i<node_chain.size(); i++)
	  {
	  cout << "\tNode chain: " << node_chain[i] << "\n";
	  //nodes[node_chain[i]].show();
	  //cout << "\n";
	  }
	*/

	//-----------------------------------------------------------------------------------------//
	//---- ACTUALLY DE-ROOT TREE (UNLIKE de_root method above that only fixes a blength=0 -----//
	//-----------------------------------------------------------------------------------------//
	// set up req'd variables
	double new_blength = -1;
	int left_of_root = nodes[the_root].get_left();
	int right_of_root = nodes[the_root].get_right();

	// If the tree has branch lengths (anything other than "-1") then
	// set the new branch between "left_of_root" and "right_of_root" equal
	// to the sum of the left and right branch lengths to the "old root".
	// This will be needed later if bool blenghts = true.
	if(nodes[left_of_root].get_blength() != -1)
	{
		new_blength = nodes[left_of_root].get_blength() + nodes[right_of_root].get_blength();
	}

	// This effectively cuts the root node out of the BTM. Can't leave it like this; the tree
	// must be re-rooted and cleaned up (code for that is below)!
	int end_of_chain = node_chain.back();               // end of node chain is the ultimate descendant in the chain
	if(end_of_chain == nodes[right_of_root].get_id())
	{
		nodes[right_of_root].set_anc(left_of_root);      // node to right of old root node has node to left of old root node as its new descendant
		nodes[right_of_root].set_blength(new_blength);   // right of old root gets the sum (or -1)
		nodes[left_of_root].set_anc(right_of_root);      // this is temp and will be updated later
	}
	else if(end_of_chain == nodes[left_of_root].get_id())
	{
		nodes[left_of_root].set_anc(right_of_root);      // node to left of old root node has node to right or old root node as its new descendant
		nodes[left_of_root].set_blength(new_blength);    // left of old root gets the sum (or -1)
		nodes[right_of_root].set_anc(left_of_root);      // this is temp and will be updated later
	}
	else
	{
		cout << "\nERROR::tree::set_root: I have fallen and I can't get up! End of chain is neither to the left or the right of the old root!\n";
		exit(0);
	}


	//----------------------------------------------------------------------------------------------//
	//---- RE-ROOT tree by resetting variables for "the_root" node relative to the new "target" ----//
	//----------------------------------------------------------------------------------------------//
	int target2 = nodes[target].get_anc();                    // new root will be BETWEEN target and anc_of_target (target2)
	nodes[target].set_anc(the_root);                          // target has a new anc (i.e., "the_root"); lf and rt descendants are unchanged

	// Now, can be one of two cases; target is to right or target is to left of its ancestor (target 2). The Operations
	// are complimentary, but one can swap right with left (free rotation on a node) without affecting anything else.
	// Hence if the case that target node is not the left of target2, I just make the swap and proceed accordingly.
	if(target == nodes[target2].get_right())                  // if case "to the right", make a swap
	{
		int new_left = nodes[target2].get_right();
		nodes[target2].set_right(nodes[target2].get_left());
		nodes[target2].set_left(new_left);
	}
	nodes[target2].set_left(nodes[target2].get_anc());        // anc of target2 copied to left descendant (old left descendant to connect to root)
	nodes[target2].set_anc(the_root);                         // "the_root" now can become new anc of target2
	nodes[the_root].set_right(target2);                       // target2 must be right descendant of root in this case
	nodes[the_root].set_left(target);                         // hence, target must be left side of root


	//------------------------------------------------------------------------------------//
	//---- now clean it up by passing over the chain of nodes that have been affected ----//
	//------------------------------------------------------------------------------------//
	for(unsigned int i=2; i<node_chain.size(); i++)                                       // Starting at node in the chain after root node
	{
		if(nodes[node_chain[i-1]].get_id() != nodes[node_chain[i]].get_anc())              // Node chain should have ancestor-descendant relationship; if not then it needs to be cleaned up
		{
			if(nodes[node_chain[i-1]].get_id() == nodes[node_chain[i]].get_left())          // Swap ancestor descendant (left in this case) relationship
			{
				int old_anc = nodes[node_chain[i]].get_anc();
				int old_lf = nodes[node_chain[i]].get_left();
				nodes[node_chain[i]].set_anc(old_lf);
				nodes[node_chain[i]].set_left(old_anc);
			}
			else if(nodes[node_chain[i-1]].get_id() == nodes[node_chain[i]].get_right())    // Swap ancestor descendant (right in this case) relationship
			{
				int old_anc = nodes[node_chain[i]].get_anc();
				int old_rt = nodes[node_chain[i]].get_right();
				nodes[node_chain[i]].set_anc(old_rt);
				nodes[node_chain[i]].set_right(old_anc);
			}
		}
	}

	/*
	  for(unsigned int i=0; i<node_chain.size(); i++)
	  {
	  cout << "Node chain (after step 2 of re-root): " << node_chain[i] << "\n";
	  nodes[node_chain[i]].show();
	  cout << "\n";
	  }
	*/

	//------------------------------------------------------------------------------------//
	//----------                 now deal with branch lengths                   ----------//
	//------------------------------------------------------------------------------------//
	// NOTE: if a node is in the "node chain", and that node is NOT connected to either the
	//       new-root nor the old-root then all the need happen is the anc b.l. be swapped
	//       with the appropriate descendant b.l.
	if(blengths) // try to correctly set blengths for the new, rooted, topology
	{
		int start = 2;
		int end = (node_chain.size()-2);
		double bl1 = nodes[node_chain[start]].get_blength();
		//cout << "bl1 = " << bl1 << "\n\n";
		double bl2;
		if(end-start >= 0)
		{
			for(int i=start; i<end; i++)                                       // Starting at node in the chain after root node
			{
				//cout << "i = " << i << "\n";
				bl2 = nodes[node_chain[i+1]].get_blength();
				//cout << "bl2 as storage = " << bl2 << "\n";
				nodes[node_chain[i+1]].set_blength(bl1);
				//cout << "bl1 before = " << bl1 << "\n";
				bl1 = bl2;
				//cout << "bl1 after = " << bl1 << "\n\n";

			}
			// Here "start" of node chain will always be node that is the right descendant of root node!
			// This node needs b.l. set to zero (violates convention req'd for b.l. estimation, but
			// those methods make the adjustment. Just in case, I also do the rotation at the end of this
			// method (if needed).
			nodes[node_chain[start]].set_blength(0);
		}
	}
	else  // set all to -1; CAUTION: This is permanent!!!
	{
		for(unsigned int i=0; i<nodes.size(); i++)
		{
			nodes[i].set_blength(-1.0);
		}
	}

	//********THIS NOW SEEMS LIKE IT IS BROKEN; CHECK IT OUT ************//

	//------------------------------------------------------------------------------------//
	//----------           now deal with the so-called "branch marks"           ----------//
	//------------------------------------------------------------------------------------//
	// NOTE: if a node is in the "node chain", and that node is NOT connected to either the
	//       new-root nor the old-root then all the need happen is the anc mark be swapped
	//       with the appropriate descendant mark.
	if(bmarks) // try to correctly set bmarks for the new, rooted, topology
	{
		int start = 2;
     	int end = (node_chain.size()-2);
     	int start_bmark = nodes[node_chain[0]].get_branch_class();
     	//cout << "Start of node chain = " << node_chain[start] << "\n";
		int bmark1 = nodes[node_chain[start]].get_branch_class();
		//cout << "bmark1 = " << bmark1 << "\n\n";
		int bmark2;
		if(end-start >= 0)
		{
			for(int i=start; i<end; i++) // Starting at node in the chain after root node
			{
				//cout << "i = " << i << "\n";
				bmark2 = nodes[node_chain[i+1]].get_branch_class();
				//cout << "bmark2 as storage = " << bmark2 << "\n";
				nodes[node_chain[i+1]].set_branch_class(bmark1);
				//cout << "bmark1 before = " << bmark1 << "\n";
				bmark1 = bmark2;
				//cout << "bmark1 after = " << bmark1 << "\n\n";
			}
			// Here "start" of node chain will always be node that is the right descendant of root node!
			// The mark at the right descendant is the mark used for the branch that has non-zero branch
			// length. To make this happen I just set both the right and left descendants = start_mark
			int rt = nodes[0].get_right();
			int lf = nodes[0].get_left();
			nodes[rt].set_branch_class(start_bmark);
			nodes[lf].set_branch_class(start_bmark);
		}
	}
	else  // set all to -1; CAUTION: This is permanent!!!
	{
		for(unsigned int i=0; i<nodes.size(); i++)
		{
			nodes[i].set_branch_class(-1.0);
		}
	}

	//------------------------------------------------------------------------------------//
	//----------              now just some final bookkeeping                   ----------//
	//------------------------------------------------------------------------------------//
	// reset all node heights to -1 before getting new node heights!!!
	for(unsigned int i=0; i<nodes.size(); i++)
	{
		nodes[i].set_height(-1);
	}

	// set the new traversal order
	//cout << "\ttree::set_root: Tree has been RE-rooted:...traversing\n";
	traverse();

	//set the new node heights
	set_node_heights();

	// Swapping if case right of root has zero b.l. in stead of left of root.
	// This maintains a convention imposed by the fast Newton-Raphson routine
	// for b.l. optimization.
	int rt = nodes[the_root].get_right();
	if(nodes[rt].get_blength() == 0.0)
	{
		int lf = nodes[the_root].get_left();
		// swap left and right
		nodes[the_root].set_right(lf);
		nodes[the_root].set_left(rt);
	}

	//all done
	/*
	  cout << "\tBTM = ";
	  euler_tour();
	  cout << "\n\n";
	*/
}


//******************************************************************************
//
//         ------------------   TreeCypher methods  --------------------
//
//******************************************************************************


//------------------------------------------------------------------------------
// public: show tree with ALL nodes numbered
//------------------------------------------------------------------------------
void tree::show_numbered_tree()
{
    mySS.str(""); //make sure to start with an empty stringstream
    
    //cout << "\t[node]:\t\t";
    int i = get_root();
    numbered_euler(i, false);
    //cout << "\n\n";
    
    string myString = mySS.str();
    cout << myString << "\n";
}


//------------------------------------------------------------------------------
// public: show tree any branch marks, or clade marks, if
//         present in tree data structure
//------------------------------------------------------------------------------
void tree::show_marked_tree()
{
    mySS.str(""); //make sure to start with an empty stringstream
    
    int i = get_root();
    marked_euler(i, false);
    
    cout << "\nNOTE: Branches in tree with labels (if any) are denoted\n";
    cout << "by a \"# mark\" in the tree file.  This marks the node to be\n";
    cout << "labelled, and the value of the label is a user-defined integer.\n\n";
    
    string myString = mySS.str();
    cout << myString << "\n";
}


//------------------------------------------------------------------------------
// public: show table of branch classes on screen (mostly for de-bugging)
//------------------------------------------------------------------------------
void tree::show_branch_class_table()
{
    cout << "\n\nBranch class details for each node in tree:\n\n";
    
    // show a node list with info relevant to setting the branch class IDs
    bool tip = false;
    cout << "\tNode  \tLeft  \tRight  \tBranch-Class \tClade-class \n";
    for(int i=0; i<num_nodes; i++)
    {
        cout << "\t" << nodes[i].get_id();  if (nodes[i].get_anc() == -1) {cout << "*";}
        if(nodes[i].is_tip()) {tip = true;}
        if(tip)
        {
            cout << "\t(Tip)";
            cout << "\t(Tip)";
        }
        else
        {
            cout << "\t" << nodes[i].get_left();
            cout << "\t" << nodes[i].get_right();
        }
        cout << "\t" << nodes[i].get_branch_class();
        cout << "\t\t" << nodes[i].get_clade_class();
        if (nodes[i].get_anc() == -1) {cout << "\t\t(*This is the root node)";}
        cout << "\n";
    }
    cout << "\n";
}


//------------------------------------------------------------------------------
// public: print tree to a file with ALL nodes numbered
//------------------------------------------------------------------------------
void tree::print_numbered_tree()
{
    mySS.str(""); //make sure to start with an empty stringstream
    
    int i = get_root();
    numbered_euler(i, true);
    
    //string myString = mySS.str();
    //cout << myString << "\n";
}


//------------------------------------------------------------------------------
// public:  print tree AND include any branch or clade marks
//------------------------------------------------------------------------------
void tree::print_marked_tree()
{
    mySS.str(""); //make sure to start with an empty stringstream
    
    int i = get_root();
    marked_euler(i, true);
}




//------------------------------------------------------------------------------
// public: returns the mark (integer) to the specified node
//------------------------------------------------------------------------------
void tree::show_node_marks(int node)
{
    cout << "\nStatus of node: " << node << "\n";
    cout << "Branch mark: " << nodes[node].get_branch_class() << "\n";
    cout << "Clade mark: " << nodes[node].get_clade_class() << "\n";
    cout << "Node label: " << nodes[node].get_label() << "\n";
    cout << "\n\n";
}



//------------------------------------------------------------------------------
// public: add the mark (integer) to the specified node
//------------------------------------------------------------------------------
void tree::set_branch_mark(int node, int mark)
{
    if(node > num_nodes || node < 0)
    {
        cout << "\n\nWARNING: the number of nodes in the tree is " << (num_nodes -1) << "!\n";
        cout << "Node 0 is used to indicate the root node. The entered value (" << node << ") is\n";
        cout << "out of bounds. This node update did NOT succeed.\n\n";
    }
    else
    {
        cout << "\n\nBranch mark UPDATED at node " << node << "...\n";
        cout << "Previous value: " << nodes[node].get_branch_class() << "\n";
        nodes[node].set_branch_class(mark);
        cout << "Current value: " << nodes[node].get_branch_class() << "\n";
    }
}



//------------------------------------------------------------------------------
// public: remove the mark (integer) to the specified node
//------------------------------------------------------------------------------
void tree::remove_branch_mark(int node)
{
    if(node > num_nodes || node < 0)
    {
        cout << "\n\nWARNING: the number of nodes in the tree is " << (num_nodes -1) << "!\n";
        cout << "Node 0 is used to indicate the root node. The entered value (" << node << ") is\n";
        cout << "out of bounds. This node update did NOT succeed.\n\n";
    }
    else
    {
        cout << "Branch mark REMOVED at node " << node << "...\n";
        cout << "Previous value: " << nodes[node].get_branch_class() << "\n";
        nodes[node].set_branch_class(0);
        cout << "Current value: " << nodes[node].get_branch_class() << "\n";
    }
}



//------------------------------------------------------------------------------
// private: The is the RECURSIVE part of the Euler tour algorithm used to get a
//          a Newick tree with each node numbered tree
//------------------------------------------------------------------------------
void tree::numbered_euler(int i, bool print)
{
    string filename = "numbered_tree.txt";
    
    //NOTE: to show a tree to the screen, simply run this method and than
    //      cout the mySS stringstream after the method done.
    
    if(nodes[i].is_tip())
    {
        //cout << nodes[i].get_label();
        //cout << "_[" << nodes[i].get_id() << "]";
        
        mySS << nodes[i].get_label();
        mySS << "_[" << nodes[i].get_id() << "]";
    }
    else
    {
        //cout << "(";                  // "on the left of node i"
        mySS << "(";
        
        numbered_euler(nodes[i].get_left(), false);   // recursion to the left
        
        //cout << ",";                  // "from below node i"
        mySS << ",";
        
        numbered_euler(nodes[i].get_right(), false);  // recursion to the right
        
        //cout << ")";                  // "on the right of node i"
        mySS << ")";
        if(i != get_root() )
        {
            //cout << nodes[i].get_id();
            mySS << nodes[i].get_id();
        }
    }
    
    if(print) // prints tree to a file
    {
        ofstream fout;
        fout.open( filename.c_str() );
        
        if(!fout.is_open())
        {
            cout << "ERROR::tree.cpp::numbered_euler: The file " << filename << "  could not be opened.\n";
            cout << "Results cannot be printed.\n";
            cout << "Something went wrong, so TreeCypher must terminate.\n";
            exit(0);
        }
        
        cout << "Filename = \"" << filename << "\"\n";
        string myString = mySS.str();
        cout << "Tree = " << myString << "\n";
        fout << myString;
        fout.close();
    }
}


//------------------------------------------------------------------------------
// private: The is the RECURSIVE part of the Euler tour algorithm used to get a
//          a Newick tree with branch marks or clade marks
//------------------------------------------------------------------------------
void tree::marked_euler(int i, bool print)
{
    int mark = 0;
    string filename = "marked_tree.txt";
    
    if(nodes[i].is_tip())
    {
        //cout << nodes[i].get_label();
        mark = nodes[i].get_branch_class();
        if(mark !=0) cout << "#" << mark;
        
        mySS << nodes[i].get_label();
        if(mark !=0) mySS << "#" << mark;
        
        
    }
    else
    {
        //cout << "(";                  // "on the left of node i"
        mySS << "(";
        
        marked_euler(nodes[i].get_left(), false);   // recursion to the left
        
        //cout << ",";                  // "from below node i"
        mySS << ",";
        
        marked_euler(nodes[i].get_right(), false);  // recursion to the right
        
        //cout << ")";                  // "on the right of node i"
        mySS << ")";
        
        if(i != get_root() )
        {
            mark = nodes[i].get_branch_class();
            if(mark !=0)
            {
                //cout << "#" << mark;
                mySS << "#" << mark;
            }
            //cout << "#" << nodes[i].get_branch_class();
        }
    }
    
    if(print)
    {
        ofstream fout;
        fout.open( filename.c_str() );
        
        if(!fout.is_open())
        {
            cout << "ERROR::tree.cpp::numbered_euler: The file " << filename << "  could not be opened.\n";
            cout << "Results cannot be printed.\n";
            cout << "Something went wrong, so TreeCypher must terminate.\n";
            exit(0);
        }
        
        cout << "Filename = \"" << filename << "\"\n";
        string myString = mySS.str();
        cout << "Tree = " << myString << "\n";
        fout << myString;
        fout.close();
    }
}


//------------------------------------------------------------------------------
// public:
//------------------------------------------------------------------------------
void tree::list_descendants(bool show, bool print)
{
    string myString;
    string filename = "monophyletic_subtrees.txt";
    mySS.str(""); //make sure to start with an empty stringstream
    
    if(show) //SHOW on screen
    {
        cout << "Internal Node   Descendant OTUs\n";
        for(int i=0; i<num_nodes; i++)
        {
            if(nodes[i].is_tip()) continue;
            cout << "          " << i << "     ";
            visit_descendants(i);
            myString = mySS.str();
            cout << myString << "\n";
            mySS.str("");
        }
    }
    if(print)
    {
        ofstream fout;
        fout.open( filename.c_str() );
        
        if(!fout.is_open())
        {
            cout << "ERROR::tree.cpp::numbered_euler: The file " << filename << "  could not be opened.\n";
            cout << "Results cannot be printed.\n";
            cout << "Something went wrong, so TreeCypher must terminate.\n";
            exit(0);
        }
        
        cout << "Filename = \"" << filename << "\"\n";
        cout << "Number of monophyletic subtrees in file = " << (num_nodes - num_otu) << "\n";
        
        fout << "Internal Node   Descendant OTUs\n";
        for(int i=0; i<num_nodes; i++)
        {
            if(nodes[i].is_tip()) continue;
            fout << "          " << i << "     ";
            visit_descendants(i);
            myString = mySS.str();
            fout << myString << "\n";
            mySS.str("");
        }
    }
}


//------------------------------------------------------------------------------
// public: the is the RECURSIVE algorithm that visits all nodes that descend
//         from a given node i; it represent at PARTIAL euler tour of a tree
//------------------------------------------------------------------------------
void tree::visit_descendants(int i)
{
    if(nodes[i].is_tip())
    {
        mySS << nodes[i].get_label() << " ";
    }
    else
    {
        visit_descendants(nodes[i].get_left());   // recursion to the left
        visit_descendants(nodes[i].get_right());  // recursion to the right
    }
}
