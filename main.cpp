/*  DendroCypher
 *  Name: main.cpp
 *  Copyright:
 *  Author: J. P. Bielawski
 *  Date: 09 Nov 2015
 *  Description: .
 */

#include <iostream>
#include <limits>
#include <fstream>
#include <math.h>
#include "matrices/matrix.h"
#include "matrices/str_matrix.h"
#include "matrices/int_matrix.h"
#include "dendrology/forestry.h"
#include "dendrology/tree.h"

using namespace std;


void show_menu();
void wait(int i);

int main()
{
    cout << "\n\n\nDendroCypher: under development.\n\n";
    cout << "\tMax number of OTUs in tree: 100\n";
    cout << "\tMax characters in OTU labels: 100\n";
    cout << "\tMultiple trees in the same file are ignored.\n";
    cout << "\tTrees must be rooted.\n";
    
    cout << "\n\n\tNOTE: This version will take as input a text file\n";
    cout << "\tcalled \"tree.txt\".  This file can contain more than\n";
    cout << "\tone tree.  The program currently works ONLY with the\n";
    cout << "\tfirst tree.  In the future the user will be able to\n";
    cout << "\tspecific which tree to work on.\n\n\n";
    
    // TreeCypher uses the tree class from Proteus.  The class treeCypher
    // is a sublcass of "tree" that implements TreeCypher specific methods.
    string filename = "tree.txt";
    tree             my_trees;
    forestry         my_forestry(filename, &my_trees);
    
    my_trees.population[0]->load_branch_marks();
    //my_trees.population[0]->load_clade_marks();
    
    int option;
    int node_id;
    int node_mark;
    string user_input;
    do {
        show_menu();
        cout << "Enter option: ";
        getline (cin,user_input);
        stringstream(user_input) >> option;
        switch (option)
        {
            case 1:
                cout << "\nShowing tree with ALL nodes numbered:\n\n";
                my_trees.population[0]->show_numbered_tree();
                wait(option);
                break;
            case 2:
                cout << "\nPrinting tree with ALL nodes numbered:\n\n";
                my_trees.population[0]->print_numbered_tree();
                wait(option);
                break;
            case 3:
                cout << "\nShowing subtrees:\n\n";
                cout << "NOTE: Because tree is rooted, subtrees are shown according\n";
                cout << "to the DESCENDANT OTUs for each internal node of the tree.\n\n";
                my_trees.population[0]->list_descendants(true, false);
                wait(option);
                break;
            case 4:
                cout << "\nPrinting subtrees according to their DESCENDANT OTUs:\n\n";
                my_trees.population[0]->list_descendants(false, true);
                wait(option);
                break;
            case 5:
                cout << "\nEnter node ID: ";
                getline (cin,user_input);
                stringstream(user_input) >> node_id;
                cout << "Enter an INTEGER to label the branch connected to this node: ";
                getline (cin,user_input);
                stringstream(user_input) >> node_mark;
                my_trees.population[0]->set_branch_mark(node_id, node_mark);
                wait(option);
                break;
            case 6:
                cout << "\nEnter node ID: ";
                getline (cin,user_input);
                stringstream(user_input) >> node_id;
                my_trees.population[0]->remove_branch_mark(node_id);
                wait(option);
                break;
            case 7:
                my_trees.population[0]->show_branch_class_table();
                wait(option);
                break;
            case 8:
                cout << "\nShowing tree with branch labels:\n\n";
                my_trees.population[0]->show_marked_tree();
                wait(option);
                break;
            case 9:
                cout << "\nPrinting tree with branch labels:\n\n";
                my_trees.population[0]->print_marked_tree();
                wait(option);
                break;
            case 10:
                break;
                
            default:
                cout << "Invalid input" << "\n";
        }
    } while (option != 10);

    my_trees.delete_trees();

   cout << "\n\nDone.\n\n";
   return 0;
}


//------------------------------------------------------------------------------
// show menu options on the screen
//------------------------------------------------------------------------------
void show_menu()
{
    cout << "\n";
    cout << "\t_______________________________________________________\n";
    cout << "\t MENU                                        Option No.\n";
    cout << "\t-------------------------------------------------------\n\n";
    cout << "\t Node IDs:\n";
    cout << "\t  + Show tree with all nodes numbered               1\n";
    cout << "\t  + Print tree to file with all nodes numbered      2\n\n";
    cout << "\t Subtrees (tree bipartitions):\n";
    cout << "\t  + Show all subtrees                               3\n";
    cout << "\t  + Print all subtrees to a file                    4\n\n";
    cout << "\t Branch and clade labels:\n";
    cout << "\t  + Add a branch label to a tree                    5\n";
    cout << "\t  + Remove a branch label from tree                 6\n";
    cout << "\t  + Show node status                                7\n";
    cout << "\t  + Show tree with branch labels                    8\n";
    cout << "\t  + Print tree with branch labels                   9\n\n";
    cout << "\t  EXIT                                             10\n";
    cout << "\t_______________________________________________________\n\n";
}

//------------------------------------------------------------------------------
// show menu options on the screen
//------------------------------------------------------------------------------
void wait(int i)
{
    string myString;
    int foobar;
    cout << "\n\n\n\tMENU: Enter the number " << i << " again to return...";
    getline (cin,myString);
    stringstream(myString) >> foobar;
}




