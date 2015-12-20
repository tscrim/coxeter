/*
  This is commands.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#include "commands.h"

#include "directories.h"
#include "error.h"
#include "fcoxgroup.h"
#include "help.h"
#include "interactive.h"
#include "special.h"
#include "typeA.h"

namespace commands {
  using namespace directories;
  using namespace error;
  using namespace fcoxgroup;
  using namespace help;
  using namespace interactive;
};

namespace {

  using namespace commands;
  using namespace stack;

  bool wgraph_warning = true;

  /* used in the definition of command trees */

  struct Empty_tag {};
  struct Interface_tag {};
  struct Main_tag {};
  struct Uneq_tag {};

  Stack<CommandTree *> treeStack;
  CoxGroup* W = 0;

  void activate(CommandTree* tree);
  void ambigAction(CommandTree* tree, const String& str);
  CommandData* ambigCommand();
  void cellCompletion(DictCell<CommandData>* cell);
  void commandCompletion(DictCell<CommandData>* cell);
  void empty_error(char* str);
  CommandTree* emptyCommandTree();
  template<class C> CommandTree* initCommandTree();
  void printCommandTree(FILE* file, DictCell<CommandData>* cell);
  void startup();

  void interface_entry();
  void interface_exit();

  void main_entry();
  void main_exit();

  void uneq_entry();
  void uneq_exit();

  void author_f();
  void betti_f();
  void coatoms_f();
  void compute_f();
  void descent_f();
  void duflo_f();
  void extremals_f();
  void fullcontext_f();
  void help_f();
  void ihbetti_f();
  void inorder_f();
  void interface_f();
  void interval_f();
  void invpol_f();
  void klbasis_f();
  void lcorder_f();
  void lcells_f();
  void lcwgraphs_f();
  void lrcorder_f();
  void lrcells_f();
  void lrcwgraphs_f();
  void lrwgraph_f();
  void lwgraph_f();
  void matrix_f();
  void mu_f();
  void not_implemented_f();
  void pol_f();
  void q_f();
  void qq_f();
  void rank_f();
  void rcorder_f();
  void rcells_f();
  void rcwgraphs_f();
  void rwgraph_f();
  void schubert_f();
  void show_f();
  void showmu_f();
  void slocus_f();
  void sstratification_f();
  void type_f();
  void uneq_f();

  const char* author_tag = "prints a message about the author";
  const char* betti_tag = "prints the ordinary betti numbers";
  const char* coatoms_tag = "prints out the coatoms of an element";
  const char* compute_tag = "prints out the normal form of an element";
  const char* descent_tag = "prints out the descent sets";
  const char* duflo_tag = "prints out the Duflo involutions";
  const char* extremals_tag =
    "prints out the k-l polynomials for the extremal pairs";
  const char* fullcontext_tag = "sets the context to the full group";
  const char* help_tag = "enters help mode";
  const char* ihbetti_tag = "prints the IH betti numbers";
  const char* input_tag = "(in help mode only) explains the input conventions";
  const char* interface_tag = "changes the interface";
  const char* interval_tag = "prints an interval in the Bruhat ordering";
  const char* intro_tag =
    "(in help mode only) prints a message for first time users";
  const char* inorder_tag = "tells whether two elements are in Bruhat order";
  const char* invpol_tag = "prints a single inverse k-l polynomial";
  const char* klbasis_tag = "prints an element of the k-l basis";
  const char* lcorder_tag = "prints the left cell order";
  const char* lcells_tag = "prints out the left k-l cells";
  const char* lcwgraphs_tag = "prints out the W-graphs of the left k-l cells";
  const char* lrcorder_tag = "prints the two-sided cell order";
  const char* lrcells_tag = "prints out the tow-sided k-l cells";
  const char* lrcwgraphs_tag =
    "prints out the W-graphs of the two-sided k-l cells";
  const char* lrwgraph_tag = "prints out the two-sided W-graph";
  const char* lwgraph_tag = "prints out the left W-graph";
  const char* matrix_tag = "prints the current Coxeter matrix";
  const char* mu_tag = "prints a single mu-coefficient";
  const char* pol_tag = "prints a single k-l polynomial";
  const char* q_tag = "exits the current mode";
  const char* qq_tag = "exits the program";
  const char* rank_tag = "resets the rank";
  const char* rcorder_tag = "prints the right cell order";
  const char* rcells_tag = "prints out the right k-l cells";
  const char* rcwgraphs_tag = "prints out the W-graphs of the right k-l cells";
  const char* rwgraph_tag = "prints out the right W-graph";
  const char* schubert_tag = "prints out the kl data for a schubert variety";
  const char* show_tag = "maps out the computation of a k-l polynomial";
  const char* showmu_tag = "maps out the computation of a mu coefficient";
  const char* slocus_tag =
    "prints the rational singular locus of the Schubert variety";
  const char* sstratification_tag =
    "prints the rational singular stratification of the Schubert variety";
  const char* type_tag =
    "resets the type and rank (hence restarts the program)";
  const char* uneq_tag = "puts the program in unequal-parameter mode";

  namespace uneq {
    void klbasis_f();
    void lcorder_f();
    void lrcorder_f();
    void lcells_f();
    void lrcells_f();
    void mu_f();
    void pol_f();
    void rcells_f();
    void rcorder_f();

    const char* lcorder_tag = "prints the left cell order";
    const char* lrcorder_tag = "prints the two-sided cell order";
    const char* lcells_tag = "prints out the left k-l cells";
    const char* lrcells_tag = "prints out the two-sided k-l cells";
    const char* mu_tag = "prints out a mu-coefficient";
    const char* pol_tag = "prints out a single k-l polynomial";
    const char* rcells_tag = "prints out the right k-l cells";
    const char* rcorder_tag = "prints the right cell order";
  };
};

namespace commands {
  void (*default_help)() = &help::default_h;
};

namespace commands {

  namespace interface {

    GroupEltInterface* in_buf = 0;

    struct In_tag {};
    struct Out_tag {};

    void in_entry();
    void in_exit();
    void out_entry();
    void out_exit();

    void abort_f();
    void alphabetic_f();
    void bourbaki_f();
    void default_f();
    void decimal_f();
    void gap_f();
    void hexadecimal_f();
    void in_f();
    void out_f();
    void permutation_f();
    void symbol_f();
    void terse_f();
    void ordering_f();

    const char* abort_tag = "leaves without modifying the interface";
    const char* alphabetic_tag = "sets alphabetic generator symbols";
    const char* bourbaki_tag = "sets Bourbaki conventions for i/o";
    const char* decimal_tag = "sets decimal generator symbols";
    const char* default_tag = "sets i/o to default mode";
    const char* gap_tag = "sets i/o to GAP mode";
    const char* hexadecimal_tag = "sets hexadecimal generator symbols";
    const char* in_tag = "enters reset-input mode";
    const char* out_tag = "enters reset-output mode";
    const char* permutation_tag =
      "sets permutation notation for i/o (in type A only)";
    const char* ordering_tag = "modifies the ordering of the generators";
    const char* terse_tag = "sets i/o to terse mode";

    namespace in {
      void alphabetic_f();
      void bourbaki_f();
      void decimal_f();
      void default_f();
      void gap_f();
      void hexadecimal_f();
      void permutation_f();
      void postfix_f();
      void prefix_f();
      void separator_f();
      void terse_f();
      const char* alphabetic_tag =
        "sets alphabetic generator symbols for input";
      const char* bourbaki_tag = "sets Bourbaki conventions for input";
      const char* decimal_tag = "sets decimal generator symbols for input";
      const char* default_tag = "sets default conventions for input";
      const char* gap_tag = "sets GAP conventions for input";
      const char* hexadecimal_tag =
        "sets hexadecimal generator symbols for input";
      const char* permutation_tag =
        "sets permutation notation for input (in type A only)";
      const char* postfix_tag = "resets the input postfix";
      const char* prefix_tag = "resets the input prefix";
      const char* separator_tag = "resets the input separator";
      const char* symbol_tag = "resets an input symbol";
      const char* terse_tag = "sets terse conventions for input";
    };

    namespace out {
      void alphabetic_f();
      void bourbaki_f();
      void decimal_f();
      void default_f();
      void gap_f();
      void hexadecimal_f();
      void permutation_f();
      void postfix_f();
      void prefix_f();
      void separator_f();
      void terse_f();
      const char* alphabetic_tag =
        "sets alphabetic generator symbols for output";
      const char* bourbaki_tag = "sets Bourbaki conventions for output";
      const char* decimal_tag = "sets decimal generator symbols for output";
      const char* default_tag = "sets default conventions for output";
      const char* gap_tag = "sets GAP conventions for output";
      const char* hexadecimal_tag =
        "sets hexadecimal generator symbols for output";
      const char* permutation_tag =
        "sets permutation notation for output (in type A only)";
      const char* postfix_tag = "resets the output postfix";
      const char* prefix_tag = "resets the output prefix";
      const char* separator_tag = "resets the output separator";
      const char* symbol_tag = "resets an output symbol";
      const char* terse_tag = "sets terse conventions for output";
    };
};

};

/*****************************************************************************

  This module contains the code for the command interface. Although overall
  I'm happy with the way it works, it suffers from a certain amount of
  clumsiness.

  The idea is that at each point in time, there is a certain active CommandTree
  object. This is basically a dictionary of recognized command names, together
  with the functions that will executed for them; in other words, something
  that should be a map in STL parlance. Actually, the active command tree is
  the top of the command tree stack treeStack; exiting the current mode means
  popping the stack; entering a new mode means pushing it on the stack.

  Each mode has an associated entry and exit function, which take care of
  initialization and clean-up duties. Actually, there is mostly one main mode;
  the entry function for this is the one which gets type and rank for the user;
  the exit function destroys the current group. Redefining type or rank means
  exiting and re-entering the main mode. In addition, there is the "empty"
  mode, active on startup only, where nothing is defined yet, and some
  auxiliary modes which temporarily hide the main mode in order to perform
  certain duties : interface mode to set the i/o preferences of the user,
  help mode for help, and also unequal-parameter mode which sets unequal
  parameters for the k-l functions; this is in fact a sort of duplicate
  main mode.

  Command completion is implemented to the extend that incomplete commands
  are recognized when non-ambiguous.

 *****************************************************************************/

/*****************************************************************************

        Chapter I -- Running the program

  This section contains the following functions :

  - ambigAction(str) : what to do with an ambiguous command;
  - mainCommandTree() : returns a pointer to the initial command tree (and
    builds it on first call);
  - relax_f() : does nothing;
  - run() : runs an interactive session;

 *****************************************************************************/

namespace commands {

void relax_f()

/*
  Does nothing.
*/

{}

void run()

/*
  This function runs an interactive session of the program.
*/

{
  static String name(0);

  activate(emptyCommandTree());

  if (ERRNO) {
    Error (ERRNO);
    return;
  }

  while (1) { /* the only way to exit from this loop is the "qq" command */
    CommandTree* tree = treeStack.top();
    tree->prompt();
    getInput(stdin,name);
    CommandData* cd = tree->find(name);
    if (cd == 0) {
      tree->error(name.ptr());
      continue;
    }
    if (cd == ambigCommand()) {
      ambigAction(tree,name);
      continue;
    }
    cd->action();
    if (cd->autorepeat) {
      tree->setAction("",cd->action);
      tree->setRepeat("",true);
    }
    else {
      tree->setAction("",&relax_f);
      tree->setRepeat("",false);
    }
  }
}

void default_error(char* str)

/*
  Default response to an unknown command.
*/

{
  Error(COMMAND_NOT_FOUND,str);
  return;
}

};

namespace {

void activate(CommandTree* tree)

/*
  Puts the tree on top of treeStack, and executes the initialization function.

*/

{
  treeStack.push(tree);
  tree->entry();

  if (ERRNO) { /* an error occured during initialization */
    Error(ERRNO);
    treeStack.pop();
    ERRNO = MODECHANGE_FAIL;
  }

  return;
}

void ambigAction(CommandTree* tree, const String& str)

/*
  Response to ambiguous commands. Prints a warning and the list of possible
  completions in the current tree on stderr.
*/

{
  static String name(0);
  bool b = true;

  print(stderr,str);
  fprintf(stderr," : ambiguous (");
  DictCell<CommandData>* cell = tree->findCell(str);
  new(&name) String(str);
  printExtensions(stderr,cell->left,name,b);
  fprintf(stderr,")\n");

  return;
}

void empty_error(char* str)

{
  CommandTree* tree = mainCommandTree();

  CommandData* cd = tree->find(str);
  if (cd == 0) {
    default_error(str);
    return;
  }
  if (cd == ambigCommand()) {
    ambigAction(tree,str);
    return;
  }
  activate(tree);
  if (ERRNO) { /* something went wrong during initialization */
    Error(ERRNO);
    return;
  }
  /* type and rank are set at this point */
  if ((cd != tree->find("type")) && (cd != tree->find("rank")))
    cd->action();
  if (cd->autorepeat) {
    tree->setAction("",cd->action);
    tree->setRepeat("",true);
  }
  else {
    tree->setAction("",&relax_f);
    tree->setRepeat("",false);
  }
  return;
}

void startup()

/*
  The response to the first carriage return. Sets the response to "help"
  to a less verbose version, and starts up the program.
*/

{
  activate(mainCommandTree());

  if (ERRNO)
    Error(ERRNO);

  return;
}

};

/*****************************************************************************

        Chapter II -- The CommandTree class.

  The purpose of a CommandTree is to get a command-name from the user (or
  perhaps from a file), and execute the corresponding command. For this,
  it maintains a tree of CommandCell s, one for each initial subword of
  each recognized command. Each CommandCell knows which command it should
  execute.

  Recognizing initial subwords allows for command completion : when the
  completion is unique, the command is executed as if the full name were
  typed. When the completion is not unique, the function for ambiguous
  commands is executed : as defined here, it prints the list of all possible
  completions in the current tree, and prompts the user again.

  The case of the empty command is special : either it does nothing, or,
  for most commands, it repeats the previous command.

  This setup supports the concept of mode : at all times, there is a current
  command tree, and in some situations this will change : new commands can
  become available, commands can change behaviour or can become unavailable.
  Help mode is an example of this. Another example is the interface command,
  which loads the interface mode tree, so that the user can set the various
  i/o parameters.

  NOTE : even though I like the actual behaviour of the setup, it is all
  rather clumsy and should be re-done. The command tree could be replaced
  with some associative container like map. The current behaviour could
  be pretty much kept as is, until we include the functionalities of
  readline.

  The following functions are defined :

   - constructors and destructors :

     - CommandTree(prompt,a,entry,error,exit,h) : builds a command tree with
       the given prompt, action a, help function h, given entry and exit
       functions, and error function (called when a command is not found);
     - ~CommandTree();

   - accessors :

     - prompt : prints the prompt;

   - manipulators :

     - add(name,tag,a,h,rep) : adds a command with the given name, tag (used
       in the online help), action a, help-action h and repetition flag rep;
     - setAction(str,a) : resets the action of the command for str;
     - setRepeat(str,b) : resets the repetition flag of the command for b;

******************************************************************************/

namespace commands {

CommandTree::CommandTree(const char* prompt,
			 void (*a)(),
			 void (*entry)(),
			 void (*error)(char*),
			 void (*exit)(),
			 void (*h)())
  :d_prompt(prompt), d_entry(entry), d_error(error), d_exit(exit)

/*
  Initializes a command tree with the given prompt and action for the
  empty command.
*/

{
  d_root->ptr = new CommandData("","",a,&relax_f,false);

  if (h) { /* add help functionality */
    d_help = new CommandTree("help",&cr_h,h);
    d_help->add("q",q_tag,&q_f,0,false);
    add("help",help_tag,&help_f,&help_h,false);
  }
}

CommandTree::~CommandTree()

/*
  The memory allocated by a CommandTree object is hidden in the dictionary
  and in the d_help pointer.
*/

{
  delete d_help;
}

/******** accessors *********************************************************/

void CommandTree::prompt() const

/*
  Prints the prompt for the command tree.
*/

{
  printf("%s : ",d_prompt.ptr());
}

/******** manipulators ******************************************************/

void CommandTree::add(const char* name, const char* tag, void (*a)(),
		      void (*h)(), bool rep)

/*
  This function adds a new command to the tree, adding new cells as
  necessary.
*/

{
  CommandData *cd = new CommandData(name,tag,a,h,rep);

  insert(name,cd);
  if (d_help && h) { /* add help functionality */
    d_help->add(name,tag,h,0,false);
  }
}

void CommandTree::setAction(const char* str, void (*a)())

/*
  Assuming that str is a fullname on the command tree, sets the response
  to str to a.

  NOTE : is a bit dangerous. Should work also when str is only a prefix.
*/

{
  CommandData* cd = find(str);
  cd->action = a;

  return;
}

void CommandTree::setRepeat(const char* str, bool b)

/*
  Assuming that str is a fullname on the command tree, sets the autorepeat
  value of the corresponding command data structure to b.
*/

{
  CommandData* cd = find(str);
  cd->autorepeat = b;

  return;
}

};

/*****************************************************************************

        Chapter III -- The CommandData class.

  The CommandData structure collects the data associated to a given command
  name. The function a defines the action associated with the command; the
  function h defines the action associated with the command in help mode.
  The flag b is set if the command should be repeated on a carriage return.

 *****************************************************************************/

namespace commands {

CommandData::CommandData(const char* const& str, const char* const& t,
			 void (*a)(), void (*h)(), bool rep)
  :name(str), tag(t), action(a), help(h), autorepeat(rep)

{}

CommandData::~CommandData()

/*
  No memory is allocated directly
*/

{}

};

/*****************************************************************************

        Chapter IV -- Building the command tree

  This section contains the functions used for the construction of the primary
  command tree, i.e., the initialization of the command module.

  The following functions are defined :

  - ambigCommand() : returns a special value flagging ambiguous commands;
  - cellCompletion(cell) : auxiliary to commandCompletion;
  - commandCompletion(tree) : finishes off the command tree;
  - emptyCommandTree() : returns a pointer to the initial command tree;
  - initCommandTree<Empty_tag> : builds the empty command tree;
  - initCommandTree<Interface_tag> : builds the interface command tree;
  - initCommandTree<Main_tag> : builds the main command tree;
  - initCommandTree<Uneq_tag> : builds the command tree for unequal parameters;
  - interfaceCommandTree() : returns a pointer to the interface command tree;
  - mainCommandTree() : returns a pointer to the main command tree;
  - uneqCommandTree() : returns a pointer to the unequal-parameter command
    tree;

 *****************************************************************************/

namespace {

CommandData* ambigCommand()

/*
  Returns a dummy command cell which is a placeholder indicating that
  ambigAction must be executed; this requires knowledge of where we are
  in the command tree.
*/

{
  static CommandData cd("","",0,0,false);
  return &cd;
}

void cellCompletion(DictCell<CommandData>* cell)

/*
  This function fills in the value fields of the cells which do not
  correspond to full names. It is assumed that the tree is traversed
  in infix (?) order, i.e. first visit left child, then cell, the right
  child, so that all longer names are already visited. This allows to
  fill in unique completions backwards, making thins easy.
*/

{
  if (cell->fullname == true)
    return;

  if (cell->uniquePrefix == false) { /* ambiguous */
    cell->ptr = ambigCommand();
    return;
  }

  if (cell->uniquePrefix == true) { /* unique completion */
    cell->ptr = cell->left->value();
    return;
  }
}

void commandCompletion(DictCell<CommandData>* cell)

/*
  This function finishes up the command tree by implementing command
  completion. This is done as follows. We traverse the command tree
  from the root. If node.fullname is true, we do nothing. Otherwise,
  if node.uniquePrefix is true, we set node.value to be equal to the
  value of the unique completion of the string recognized by node in
  the dictionary. If node.uniquePrefix is false, we set node.value
  to ambigCommand().
*/

{
  if (cell == 0)
    return;

  commandCompletion(cell->left);
  cellCompletion(cell);
  commandCompletion(cell->right);
}

template<> CommandTree* initCommandTree<Empty_tag>()

/*
  This function builds the initial command tree of the program. The idea
  is that all commands on the main command tree will be considered entry
  commands, and so will do the necessary initialization. This is achieved
  thru the special error function.
*/

{
  static CommandTree tree("coxeter",&startup,&relax_f,&empty_error,&relax_f,
			  &intro_h);

  tree.add("author","author_tag",&author_f,&relax_f,false);
  tree.add("qq",qq_tag,&qq_f,&qq_h,false);

  commandCompletion(tree.root());

  tree.helpMode()->add("intro",intro_tag,&intro_h,0,false);

  commandCompletion(tree.helpMode()->root());

  return &tree;
}

CommandTree* emptyCommandTree()

/*
  Returns a pointer to the initial command tree of the program, building it on
  the first call.
*/

{
  static CommandTree* tree = initCommandTree<Empty_tag>();
  return tree;
}

template<> CommandTree* initCommandTree<Interface_tag>()

/*
  This function builds the interface command tree; this makes available the
  various little commands that are needed to reset the interface, and that
  have no reason to be clogging up the main command tree.
*/

{
  static CommandTree tree("interface",&relax_f,&interface_entry,&default_error,
			  &interface_exit,&interface_help);

  tree.add("alphabetic",commands::interface::alphabetic_tag,
	   &commands::interface::alphabetic_f,&help::interface::alphabetic_h);
  tree.add("bourbaki",commands::interface::bourbaki_tag,
	   &commands::interface::bourbaki_f,&help::interface::bourbaki_h);
  tree.add("decimal",commands::interface::decimal_tag,
	   &commands::interface::decimal_f,&help::interface::decimal_h);
  tree.add("default",commands::interface::default_tag,
	   &commands::interface::default_f,&help::interface::default_h);
  tree.add("gap",commands::interface::out::gap_tag,
	   &commands::interface::out::gap_f, &help::interface::gap_h);
  tree.add("hexadecimal",commands::interface::hexadecimal_tag,
	   &commands::interface::hexadecimal_f,
	   &help::interface::hexadecimal_h);
  tree.add("in",commands::interface::in_tag,&commands::interface::in_f,
	   help::interface::in_h,false);
  tree.add("ordering",commands::interface::ordering_tag,
	   &commands::interface::ordering_f,help::interface::ordering_h,false);
  tree.add("out",commands::interface::out_tag,&commands::interface::out_f,
	   help::interface::out_h,false);
  tree.add("permutation",commands::interface::permutation_tag,
	   &commands::interface::permutation_f,
	   &help::interface::permutation_h);
  tree.add("q",q_tag,&q_f,0,false);
  tree.add("terse",commands::interface::out::terse_tag,
	   &commands::interface::out::terse_f, &help::interface::out::terse_h);

  commandCompletion(tree.root());
  commandCompletion(tree.helpMode()->root());

  return &tree;
}

template<> CommandTree* initCommandTree<commands::interface::In_tag>()

/*
  This function builds the command tree for the input-modification mode.
*/

{
  using namespace commands::interface;

  static CommandTree tree("in",&relax_f,&in_entry,&default_error,
			  &in_exit,&help::interface::in_help);

  tree.add("q",q_tag,&q_f,0,false);

  tree.add("abort",abort_tag,&abort_f,&help::interface::abort_h);
  tree.add("alphabetic",in::alphabetic_tag,&in::alphabetic_f,
	   &help::interface::in::alphabetic_h,false);
  tree.add("bourbaki",in::bourbaki_tag,&in::bourbaki_f,
	   &help::interface::in::bourbaki_h);
  tree.add("decimal",in::decimal_tag,&in::decimal_f,
	   &help::interface::in::decimal_h,false);
  tree.add("default",in::default_tag,&in::default_f,
	   &help::interface::in::default_h);
  tree.add("gap",in::gap_tag,&in::gap_f,&help::interface::in::gap_h);
  tree.add("hexadecimal",in::hexadecimal_tag,&in::hexadecimal_f,
	   &help::interface::in::hexadecimal_h,false);
  tree.add("permutation",in::permutation_tag,&in::permutation_f,
	   &help::interface::in::permutation_h,false);
  tree.add("postfix",in::postfix_tag,&in::postfix_f,
	   &help::interface::in::postfix_h);
  tree.add("prefix",in::prefix_tag,&in::prefix_f,
	   &help::interface::in::prefix_h);
  tree.add("separator",in::separator_tag,
	   &in::separator_f,&help::interface::in::separator_h);
  tree.add("symbol",in::symbol_tag,&symbol_f,
	   &help::interface::in::symbol_h);
  tree.add("terse",in::terse_tag,&in::terse_f,&help::interface::in::terse_h);

  commandCompletion(tree.root());
  commandCompletion(tree.helpMode()->root());

  return &tree;
}

template<> CommandTree* initCommandTree<commands::interface::Out_tag>()

/*
  This function builds the command tree for the output-modification mode.
*/

{
  using namespace commands::interface;

  static CommandTree tree("out",&relax_f,&out_entry,&default_error,
			  &out_exit,&help::interface::out_help);

  tree.add("q",q_tag,&q_f,0,false);

  tree.add("alphabetic",out::alphabetic_tag,&out::alphabetic_f,
	   &help::interface::out::alphabetic_h,false);
  tree.add("bourbaki",out::bourbaki_tag,&out::bourbaki_f,
	   &help::interface::out::bourbaki_h);
  tree.add("decimal",out::decimal_tag,&out::decimal_f,
	   &help::interface::out::decimal_h,false);
  tree.add("default",out::default_tag,&out::default_f,
	   &help::interface::out::default_h);
  tree.add("gap",out::gap_tag,&out::gap_f,
	   &help::interface::out::gap_h);
  tree.add("hexadecimal",out::hexadecimal_tag,&out::hexadecimal_f,
	   &help::interface::out::hexadecimal_h,false);
  tree.add("permutation",out::permutation_tag,&out::permutation_f,
	   &help::interface::out::permutation_h,false);
  tree.add("postfix",out::postfix_tag,&out::postfix_f,
	   &help::interface::out::postfix_h);
  tree.add("prefix",out::prefix_tag,&out::prefix_f,
	   &help::interface::out::prefix_h);
  tree.add("separator",out::separator_tag,
	   &out::separator_f,&help::interface::out::separator_h);
  tree.add("symbol",out::symbol_tag,&symbol_f,
	   &help::interface::out::symbol_h);
  tree.add("terse",out::terse_tag,&out::terse_f,
	   &help::interface::out::terse_h);

  commandCompletion(tree.root());
  commandCompletion(tree.helpMode()->root());

  return &tree;
}

};

namespace commands {

CommandTree* interface::inCommandTree()

{
  static CommandTree* tree = initCommandTree<In_tag>();
  return tree;
}

CommandTree* interface::outCommandTree()

{
  static CommandTree* tree = initCommandTree<Out_tag>();
  return tree;
}

CommandTree* interfaceCommandTree()

/*
  Returns a pointer to the interface command tree, building it on the first
  call.
*/

{
  static CommandTree* tree = initCommandTree<Interface_tag>();
  return tree;
}

};

namespace {

template<> CommandTree* initCommandTree<Main_tag>()

/*
  This function builds the main command tree, the one that is being run on
  startup. Auxiliary trees may be grafted onto this one (thru the pushdown
  stack treeStack) by some functions needing to be in special modes.
*/

{
  static CommandTree tree("coxeter",&relax_f,&main_entry,&default_error,
			  &main_exit,&main_help);

  tree.add("author",author_tag,&author_f,&relax_f,false);
  tree.add("betti",betti_tag,&betti_f,&betti_h,false);
  tree.add("coatoms",coatoms_tag,&coatoms_f,&coatoms_h);
  tree.add("compute",compute_tag,&compute_f,&compute_h);
  tree.add("descent",descent_tag,&descent_f,&descent_h);
  tree.add("duflo",duflo_tag,&duflo_f,&duflo_h);
  tree.add("extremals",extremals_tag,&extremals_f,&extremals_h);
  tree.add("fullcontext",fullcontext_tag,&fullcontext_f,&fullcontext_h);
  tree.add("ihbetti",ihbetti_tag,&ihbetti_f,&ihbetti_h,false);
  tree.add("interface",interface_tag,&interface_f,&interface_h,false);
  tree.add("interval",interval_tag,&interval_f,&interval_h,false);
  tree.add("inorder",inorder_tag,&inorder_f,&inorder_h);
  tree.add("invpol",invpol_tag,&invpol_f,&invpol_h);
  tree.add("lcorder",lcorder_tag,&lcorder_f,&lcorder_h,false);
  tree.add("lcells",lcells_tag,&lcells_f,&lcells_h,false);
  tree.add("lcwgraphs",lcwgraphs_tag,&lcwgraphs_f,&lcwgraphs_h,false);
  tree.add("lrcorder",lrcorder_tag,&lrcorder_f,&lrcorder_h,false);
  tree.add("lrcells",lrcells_tag,&lrcells_f,&lrcells_h,false);
  tree.add("lrcwgraphs",lrcwgraphs_tag,&lrcwgraphs_f,&lrcwgraphs_h,false);
  tree.add("lrwgraph",lrwgraph_tag,&lrwgraph_f,&lrwgraph_h,false);
  tree.add("lwgraph",lwgraph_tag,&lwgraph_f,&lwgraph_h,false);
  tree.add("klbasis",klbasis_tag,&klbasis_f,&klbasis_h,true);
  tree.add("matrix",matrix_tag,&matrix_f,&matrix_h);
  tree.add("mu",mu_tag,&mu_f,&mu_h);
  tree.add("pol",pol_tag,&pol_f,&pol_h);
  tree.add("q",q_tag,&q_f,0,false);
  tree.add("qq",qq_tag,&qq_f,&qq_h,false);
  tree.add("rank",rank_tag,&rank_f,&rank_h,false);
  tree.add("rcorder",rcorder_tag,&rcorder_f,&rcorder_h,false);
  tree.add("rcells",rcells_tag,&rcells_f,&rcells_h,false);
  tree.add("rcwgraphs",rcwgraphs_tag,&rcwgraphs_f,&rcwgraphs_h,false);
  tree.add("rwgraph",rwgraph_tag,&rwgraph_f,&rwgraph_h,false);
  tree.add("schubert",schubert_tag,&schubert_f,&schubert_h);
  tree.add("show",show_tag,&show_f,&show_h);
  tree.add("showmu",showmu_tag,&showmu_f,&showmu_h);
  tree.add("slocus",slocus_tag,&slocus_f,&slocus_h);
  tree.add("sstratification",sstratification_tag,&sstratification_f,
	   &sstratification_h);
  tree.add("type",type_tag,&type_f,&type_h,false);
  tree.add("uneq",uneq_tag,&uneq_f,&uneq_h,false);

  special::addSpecialCommands(&tree);

  commandCompletion(tree.root());

  tree.helpMode()->add("intro",intro_tag,&intro_h,0,false);
  tree.helpMode()->add("input",input_tag,&input_h,0,false);

  commandCompletion(tree.helpMode()->root());

  return &tree;
}

};

namespace commands {

CommandTree* mainCommandTree()

/*
  Returns a pointer to the main command tree of the program, building it on
  the first call.
*/

{
  static CommandTree* tree = initCommandTree<Main_tag>();
  return tree;
}

};

namespace {

template<> CommandTree* initCommandTree<Uneq_tag>()

/*
  This function builds the unequal-parameter command tree. It contains
  essentially the same functions as the main command tree, except that the
  unequal-parameter versions have been substituted for the k-l functions.
*/

{
  static CommandTree tree("uneq",&relax_f,&uneq_entry,&default_error,
			  &uneq_exit,&uneq_help);

  tree.add("author",author_tag,&author_f,&relax_f,false);
  tree.add("coatoms",coatoms_tag,&coatoms_f,&coatoms_h);
  tree.add("compute",compute_tag,&compute_f,&compute_h);
  tree.add("descent",descent_tag,&descent_f,&descent_h);
  tree.add("fullcontext",fullcontext_tag,&fullcontext_f,&fullcontext_h);
  tree.add("interface",interface_tag,&interface_f,&interface_h,false);
  tree.add("klbasis",klbasis_tag,&uneq::klbasis_f,&help::uneq::klbasis_h,true);
  tree.add("lcorder",uneq::lcorder_tag,&uneq::lcorder_f,
	   &help::uneq::lcorder_h,false);
  tree.add("lrcorder",uneq::lrcorder_tag,&uneq::lrcorder_f,
	   &help::uneq::lrcorder_h,false);
  tree.add("lcells",uneq::lcells_tag,&uneq::lcells_f,&help::uneq::lcells_h,
	   false);
  tree.add("lrcells",uneq::lrcells_tag,&uneq::lrcells_f,&help::uneq::lrcells_h,
	   false);
  tree.add("matrix",matrix_tag,&matrix_f,&matrix_h);
  tree.add("mu",uneq::mu_tag,&uneq::mu_f,&help::uneq::mu_h);
  tree.add("pol",uneq::pol_tag,&uneq::pol_f,&help::uneq::pol_h);
  tree.add("rcells",uneq::rcells_tag,&uneq::rcells_f,&help::uneq::rcells_h,
	   false);
  tree.add("rcorder",uneq::rcorder_tag,&uneq::rcorder_f,
	   &help::uneq::rcorder_h,false);
  tree.add("q",q_tag,&q_f,0,false);
  tree.add("qq",qq_tag,&qq_f,&qq_h,false);

  commandCompletion(tree.root());
  commandCompletion(tree.helpMode()->root());

  return &tree;
}

};

namespace commands {

CommandTree* uneqCommandTree()

/*
  Returns a pointer to the uneq command tree, building it on the first call.
*/

{
  static CommandTree* tree = initCommandTree<Uneq_tag>();
  return tree;
}

};

/*****************************************************************************

        Chapter V -- Functions for the predefined commands.

  This section contains the functions defining the responses to the various
  commands which are provided by the program. The functions are placed in
  the unnamed namespace defined in this file.

  The following functions are defined :

  - author_f() : response to "author";
  - betti_f() : response to "betti";
  - coatoms_f() : response to "coatoms";
  - compute_f() : response to "compute";
  - descent_f() : response to "descent";
  - duflo_f() : response to "duflo";
  - extremals_f() : response to "extremals";
  - fullcontext_f() : response to "fullcontext";
  - ihbetti_f() : response to "ihbetti";
  - inorder_f() : response to "inorder";
  - interface_f() : response to "interface";
  - help_f() : repsonse to "help";
  - klbasis_f() : response to "klbasis";
  - lcorder_f() : response to "lcorder";
  - lcells_f() : response to "lcells";
  - lcwgraphs_f() : response to "lcwgraphs";
  - lrcorder_f() : response to "lrcorder";
  - lrcells_f() : response to "lrcells";
  - lrcwgraphs_f() : response to "lrcwgraphs";
  - lrwgraph_f() : response to "lrwgraph";
  - lwgraph_f() : response to "lwgraph";
  - matrix_f() : response to "matrix";
  - mu_f() : response to "mu";
  - not_implemented_f() : response for not (yet) implemented features;
  - q_f() : response to "q";
  - qq_f() : response to "qq";
  - rank_f() : response to "rank";
  - rcorder_f() : response to "rcorder";
  - rcells_f() : response to "rcells";
  - rcwgraphs_f() : response to "rcwgraphs";
  - rwgraph_f() : response to "rwgraph";
  - relax_f() : does nothing;
  - schubert_f() : response to "schubert";
  - show_f() : response to "show";
  - showmu_f() : response to "showmu";
  - slocus_f() : response to "slocus";
  - sstratification_f() : response to "sstratification";
  - type_f() : response to "type";
  - uneq_f() : response to "uneq";

  In uneq mode :

  - klbasis_f() : response to "klbasis";
  - lcorder_f() : response to "lcorder";
  - lcells_f() : response to "lcells";
  - lrcorder_f() : response to "lrcorder";
  - lrcells_f() : response to "lrcells";
  - mu_f() : response to "mu";
  - pol_f() : response to "pol";
  - rcells_f() : response to "rcells";
  - rcorder_f() : response to "rcorder";

  In interface mode :

  - abort_f() : aborts input interface modification;
  - alphabetic_f() : sets alphabetic generator symbols;
  - bourbaki_f() : sets bourbaki conventions;
  - decimal_f() : sets decimal generator symbols;
  - default_f() : sets default i/o;
  - gap_f() : sets GAP-style i/o;
  - hexadecimal_f() : sets hexadecimal generator symbols;
  - ordering_f() : changes the ordering of the generators;
  - symbol_f() : resets generator symbol;
  - terse_f() : sets terse style i/o;
  - in::alphabetic_f() : sets alphabetic generator symbols;
  - in::bourbaki_f() : sets bourbaki conventions;
  - in::decimal_f() : sets decimal generator symbols;
  - in::default_f() : sets default-style input;
  - in::gap_f() : sets GAP-style input;
  - in::hexadecimal_f() : sets hexadecimal generator symbols;
  - in::permutation_f() : sets permutation input;
  - in::postfix_f() : resets input postfix;
  - in::prefix_f() : resets input prefix;
  - in::separator_f() : resets input separator;
  - in::terse_f() : sets terse-style input;
  - out::alphabetic_f() : sets alphabetic generator symbols;
  - out::bourbaki_f() : sets bourbaki conventions;
  - out::decimal_f() : sets decimal generator symbols;
  - out::default_f() : sets default-style output;
  - out::gap_f() : sets GAP-style output;
  - out::hexadecimal_f() : sets hexadecimal generator symbols;
  - out::permutation_f() : sets permutation output;
  - out::postfix_f() : resets output postfix;
  - out::prefix_f() : resets output prefix;
  - out::separator_f() : resest output separator;
  - out::terse_f() : sets terse-style output;

 *****************************************************************************/

namespace {

void author_f()

/*
  Prints a message about the author.
*/

{
  printFile(stderr,"author.mess",MESSAGE_DIR);
  return;
}

void betti_f()

/*
  Prints out the ordinary betti numbers of [e,y].

  NOTE : could be *much* improved! In particular, we would want to have
  betti(x,y) for x <= y.
*/

{
  static CoxWord g(0);

  printf("enter your element (finish with a carriage return) :\n");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputTraits& traits = W->outputTraits();
  printBetti(stdout,y,W->schubert(),traits);

  return;
}

void coatoms_f()

/*
  Prints out the coatoms of a given element, computing them in elementary
  fashion.
*/

{
  static CoxWord g(0);

  printf("enter your element (finish with a carriage return) :\n");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  List<CoxWord> c(0);
  W->coatoms(c,g);

  for (Ulong j = 0; j < c.size(); ++j) {
    W->print(stdout,c[j]);
    printf("\n");
  }

  return;
}

void compute_f()

/*
  Gets an element from the user, and prints out its normal form.
*/

{
  static CoxWord g(0);

  printf("enter your element (finish with a carriage return) :\n");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  W->normalForm(g);
  W->print(stdout,g);
  if (SmallCoxGroup* Ws = dynamic_cast<SmallCoxGroup*> (W)) {
    CoxNbr x = 0;
    Ws->prodD(x,g);
    printf(" (#%lu)",static_cast<Ulong>(x));
  }
  CoxNbr x = W->contextNumber(g);
  if (x != undef_coxnbr)
    printf(" (%s%lu)","%",static_cast<Ulong>(x));
  printf("\n");

  return;
}

void descent_f()

/*
  Prints the left and right descent sets.
*/

{
  static CoxWord g(0);

  printf("enter your element (finish with a carriage return) :\n");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  LFlags f = W->ldescent(g);
  printf("L:");
  W->printFlags(stdout,f);
  printf("; R:");
  f = W->rdescent(g);
  W->printFlags(stdout,f);
  printf("\n");

  return;
}

void duflo_f()

/*
  Prints the Duflo involutions. Works for finite groups only.
*/

{
  if (!isFiniteType(W)) {
    printFile(stderr,"duflo.mess",MESSAGE_DIR);
    return;
  }

  FiniteCoxGroup* Wf = dynamic_cast<FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputFile file;
  OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),dufloH,traits);
  printDuflo(file.f(),Wf->duflo(),Wf->lCell(),Wf->kl(),W->interface(),traits);

  return;
}

void extremals_f()

/*
  Prints out the list of extremal elements x <= y (this is part of the schubert
  command.
*/

{
  static CoxWord g(0);

  printf("Enter your element (finish with a carriage-return) :\n");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),extremalsH,traits);
  printExtremals(file.f(),y,W->kl(),W->interface(),traits);

  return;
}

void fullcontext_f()

/*
  Response to the fullcontext command. This sets the context to the full
  group. Of course, it works only for finite groups.
*/

{
  if (!isFiniteType(W)) {
    printFile(stderr,"fullcontext.mess",MESSAGE_DIR);
    return;
  }

  FiniteCoxGroup* Wf = dynamic_cast<FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
  }

  return;
}

void help_f()

/*
  Response to the help command.
*/

{
  activate(treeStack.top()->helpMode());
  return;
}

void ihbetti_f()

/*
  Prints out the IH betti numbers of [e,y].
*/

{
  static CoxWord g(0);

  printf("enter your element (finish with a carriage return) :\n");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputTraits& traits = W->outputTraits();
  printIHBetti(stdout,y,W->kl(),traits);

  return;
}

void interface_f()

/*
  Response to the interface command.
*/

{
  activate(interfaceCommandTree());
  return;
}

void interval_f()

/*
  Response to the interval command.
*/

{
  CoxWord g(0);
  CoxWord h(0);

  fprintf(stdout,"first : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  fprintf(stdout,"second : ");
  h = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (not W->inOrder(g,h)) {
    fprintf(stderr,"the two elements are not in order\n");
    return;
  }

  W->extendContext(h);

  CoxNbr x = W->contextNumber(g);
  CoxNbr y = W->contextNumber(h);

  OutputFile file;

  BitMap b(W->contextSize());
  W->extractClosure(b,y);

  BitMap::ReverseIterator b_rend = b.rend();
  List<CoxNbr> res(0);

  for (BitMap::ReverseIterator i = b.rbegin(); i != b_rend; ++i)
    if (not W->inOrder(x,*i)) {
      BitMap bi(W->contextSize());
      W->extractClosure(bi,*i);
      CoxNbr z = *i; // andnot will invalidate iterator
      b.andnot(bi);
      b.setBit(z);   // otherwise the decrement will not be correct
    } else
      res.append(*i);

  schubert::NFCompare nfc(W->schubert(),W->ordering());
  Permutation a(res.size());
  sortI(res,nfc,a);

  for (size_t j = 0; j < res.size(); ++j) {
    W->print(file.f(),res[a[j]]);
    fprintf(file.f(),"\n");
  }

  return;
}

/**
  Response to the inorder command. This will tell whether two elements
  are comparable in Bruhat order, using only the elementary string operations
  (and hence not consuming any memory.)
*/
void inorder_f()
{
  CoxWord g(0);
  CoxWord h(0);
  List<Length> a(0);

  fprintf(stdout,"first : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  fprintf(stdout,"second : ");
  h = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (W->inOrder(a,g,h)) {
    fprintf(stdout,"true :   ");
    Ulong i = 0;
    for (Ulong j = 0; j < a.size(); ++j) {
      while (i < a[j]) {
	W->printSymbol(stdout,h[i]-1); // Conversion CoxLetter -> Generator
	++i;
      }
      fprintf(stdout,".");
      ++i;
    }
    while (i < h.length()) {
      W->printSymbol(stdout,h[i]-1);
      ++i;
    }
    fprintf(stdout,"\n");
  }
  else
    fprintf(stdout,"false\n");
}

void invpol_f()

/*
  Response to the invpol command. This prints out a single inverse k-l
  polynomial, without details.
*/

{
  CoxWord g(0);

  fprintf(stdout,"first : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  CoxNbr x = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  fprintf(stdout,"second : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (!W->inOrder(x,y)) {
    fprintf(stderr,"the two elements are not in Bruhat order\n");
    return;
  }

  const invkl::KLPol& pol = W->invklPol(x,y);
  if (ERRNO) {
    Error(ERRNO,x,y);
    return;
  }

  print(stdout,pol,"q");
  printf("\n");

  return;
}

void klbasis_f()

/*
  Prints out one element in the k-l basis of the group, in the format
  defined by the current output mode.
*/

{
  CoxWord g(0);

  printf("enter your element (finish with a carriage return) :\n");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  kl::HeckeElt h(0);

  W->cBasis(h,y);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),basisH,traits);
  printAsBasisElt(file.f(),h,W->schubert(),W->interface(),traits);

  return;
}

void lcorder_f()

/*
  Prints the left cell order of the current group in a file. Works only for
  finite groups.
*/

{
  if (!isFiniteType(W)) {
    printFile(stderr,"lcorder.mess",MESSAGE_DIR);
    return;
  }

  FiniteCoxGroup* Wf = dynamic_cast<FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputFile file;
  OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),lCOrderH,traits);
  printLCOrder(file.f(),Wf->kl(),Wf->interface(),traits);

  return;
}

void lcells_f()

/*
  This function prints out the left cells in the group.
*/

{
  if (!isFiniteType(W)) {
    printFile(stderr,"lcells.mess",MESSAGE_DIR);
    return;
  }

  FiniteCoxGroup* Wf = dynamic_cast<FiniteCoxGroup*> (W);

  OutputFile file;
  OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),lCellsH,traits);
  printLCells(file.f(),Wf->lCell(),Wf->kl(),Wf->interface(),traits);

  return;
}

void lcwgraphs_f()

/*
  This function prints out the W-graphs of the left cells in the group.
  It works only for finite groups currently.
*/

{
  if (!isFiniteType(W)) {
    printFile(stderr,"lcells.mess",MESSAGE_DIR);
    return;
  }

  FiniteCoxGroup* Wf = dynamic_cast<FiniteCoxGroup*> (W);

  OutputFile file;
  OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),lCellWGraphsH,traits);
  printLCellWGraphs(file.f(),Wf->lCell(),Wf->kl(),W->interface(),traits);

  return;
}

void lrcorder_f()

/*
  Prints the two-sided cell order of the current group in a file. Works only
  for finite groups.
*/

{
  if (!isFiniteType(W)) {
    printFile(stderr,"lrcorder.mess",MESSAGE_DIR);
    return;
  }

  FiniteCoxGroup* Wf = dynamic_cast<FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputFile file;
  OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),lrCOrderH,traits);
  printLRCOrder(file.f(),Wf->kl(),Wf->interface(),traits);

  return;
}

void lrcells_f()

/*
  This function prints out the two-sided cells in the group, together with the
  corresponding W-graphs. Works only for finite groups.
*/

{
  if (!isFiniteType(W)) {
    printFile(stderr,"lrcells.mess",MESSAGE_DIR);
    return;
  }

  FiniteCoxGroup* Wf = dynamic_cast<FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputFile file;
  OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),lrCellsH,traits);
  printLRCells(file.f(),Wf->lrCell(),Wf->kl(),Wf->interface(),traits);

  return;
}

void lrcwgraphs_f()

/*
  This function prints out the W-graphs of the two-sided cells in the group.
  It works only for finite groups currently.
*/

{
  if (!isFiniteType(W)) {
    printFile(stderr,"lcells.mess",MESSAGE_DIR);
    return;
  }

  FiniteCoxGroup* Wf = dynamic_cast<FiniteCoxGroup*> (W);

  OutputFile file;
  OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),lrCellWGraphsH,traits);
  printLRCellWGraphs(file.f(),Wf->lrCell(),Wf->kl(),W->interface(),traits);

  return;
}

void lrwgraph_f()

/*
  Prints out the two-sided wgraph of the current context.
*/

{
  if (!W->isFullContext() && wgraph_warning) {
    printFile(stderr,"wgraph.mess",MESSAGE_DIR);
    printf("continue ? y/n\n");
    if (!yesNo())
      return;
    printf("print this message next time ? y/n\n");
    if (!yesNo())
      wgraph_warning = false;
  }

  W->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputFile file;
  OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),lrWGraphH,traits);
  printLRWGraph(file.f(),W->kl(),W->interface(),traits);

  return;
}

void lwgraph_f()

/*
  Prints out the left wgraph of the current context.
*/

{
  if (!W->isFullContext() && wgraph_warning) {
    printFile(stderr,"wgraph.mess",MESSAGE_DIR);
    printf("continue ? y/n\n");
    if (!yesNo())
      return;
    printf("print this message next time ? y/n\n");
    if (!yesNo())
      wgraph_warning = false;
  }

  W->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputFile file;
  OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),lWGraphH,traits);
  printLWGraph(file.f(),W->kl(),W->interface(),traits);

  return;
}

void matrix_f()

/*
  Prints the Coxeter matrix.
*/

{
  interactive::printMatrix(stdout,W);

  return;
}

void not_implemented_f()

/*
  Response for not (yet) implemented commands.
*/

{
  fprintf(stderr,"Sorry, not implemented yet\n");
  return;

}

void mu_f()

/*
  Response to the mu command. This prints out a single mu-coefficient,
  without details.
*/

{
  static CoxWord g(0);

  fprintf(stdout,"first : ");
  g = interactive::getCoxWord(W);
  CoxNbr x = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  fprintf(stdout,"second : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (!W->inOrder(x,y)) {
    fprintf(stderr,"the two elements are not in Bruhat order\n");
    return;
  }

  KLCoeff mu = W->mu(x,y);
  if (ERRNO) {
    Error(ERRNO,x,y);
    return;
  }

  printf("%lu\n",static_cast<Ulong>(mu));

  return;
}

void pol_f()

/*
  Response to the pol command. This prints out a single polynomial, without
  details.
*/

{
  static CoxWord g(0);

  fprintf(stdout,"first : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  CoxNbr x = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  fprintf(stdout,"second : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (!W->inOrder(x,y)) {
    fprintf(stderr,"the two elements are not in Bruhat order\n");
    return;
  }

  const kl::KLPol& pol = W->klPol(x,y);
  if (ERRNO) {
    Error(ERRNO,x,y);
    return;
  }

  print(stdout,pol,"q");
  printf("\n");

  return;
}

void q_f()

/*
  Exits the current mode. If there is a problem on exit, the exit function
  has the option of setting an error, thus preventing the exit.
*/

{
  CommandTree* tree = treeStack.top();
  tree->exit();

  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  treeStack.pop();

  return;
}

void qq_f()

/*
  Exits the program.
*/

{
  while(treeStack.size()) {
    CommandTree* tree = treeStack.top();
    tree->exit();
    treeStack.pop();
  }

  exit(0);
}

void rank_f()

/*
  Sets the rank.
*/

{
  CoxGroup* Wloc = interactive::allocCoxGroup(W->type());

  if (ERRNO) {
    Error(ERRNO);
  }
  else {
    W = Wloc;
  }

  return;
}

void rcorder_f()

/*
  Prints the right cell order of the current group in a file. Works only for
  finite groups.
*/

{
  if (!isFiniteType(W)) {
    printFile(stderr,"rcorder.mess",MESSAGE_DIR);
    return;
  }

  FiniteCoxGroup* Wf = dynamic_cast<FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputFile file;
  OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),rCOrderH,traits);
  printRCOrder(file.f(),Wf->kl(),Wf->interface(),traits);

  return;
}

void rcells_f()

/*
  This function prints out the right cells in the group, together with the
  corresponding W-graphs.
*/

{
  if (!isFiniteType(W)) {
    printFile(stderr,"rcells.mess",MESSAGE_DIR);
    return;
  }

  FiniteCoxGroup* Wf = dynamic_cast<FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputFile file;
  OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),rCellsH,traits);
  printRCells(file.f(),Wf->rCell(),Wf->kl(),Wf->interface(),traits);

  return;
}

void rcwgraphs_f()

/*
  This function prints out the W-graphs of the right cells in the group.
  It works only for finite groups currently.
*/

{
  if (!isFiniteType(W)) {
    printFile(stderr,"lcells.mess",MESSAGE_DIR);
    return;
  }

  FiniteCoxGroup* Wf = dynamic_cast<FiniteCoxGroup*> (W);

  OutputFile file;
  OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),rCellWGraphsH,traits);
  printRCellWGraphs(file.f(),Wf->rCell(),Wf->kl(),W->interface(),traits);

  return;
}

void rwgraph_f()

/*
  Prints out the right wgraph of the current context.
*/

{
  if (!W->isFullContext() && wgraph_warning) {
    printFile(stderr,"wgraph.mess",MESSAGE_DIR);
    printf("continue ? y/n\n");
    if (!yesNo())
      return;
    printf("print this message next time ? y/n\n");
    if (!yesNo())
      wgraph_warning = false;
  }

  W->fillMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputFile file;
  OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),rWGraphH,traits);
  printRWGraph(file.f(),W->kl(),W->interface(),traits);

  return;
}

void schubert_f()

/*
  Response to the schubert command. This will print out the information
  corresponding to one element in the k-l basis, and the information on
  the singularities of the corresponding Schubert variety, in the format
  popularized by Goresky.
*/

{
  static CoxWord g(0);

  printf("Enter your element (finish with a carriage-return) :\n");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),closureH,traits);
  printClosure(file.f(),y,W->kl(),W->interface(),traits);

  return;
}

void show_f()

/*
  Response to the show command. This maps out the computation of a Kazhdan-
  Lusztig polynomial, letting you choose the descent generator.

  If no generator is given, the default descent generator is used.
*/

{
  static CoxWord g(0);

  fprintf(stdout,"first : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  CoxNbr x = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  fprintf(stdout,"second : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (!W->inOrder(x,y)) {
    fprintf(stderr,"the two elements are not in Bruhat order\n");
    return;
  }

  fprintf(stdout,"generator (carriage return for default) : ");
  LFlags f = W->descent(y);
  Generator s = interactive::getGenerator(W,f);
  if (ERRNO) {
    Error (ERRNO);
    return;
  }

  interactive::OutputFile file;
  showKLPol(file.f(),W->kl(),x,y,W->interface(),s);

  return;
}

void showmu_f()

/*
  Response to the showmu command. This maps out the computation of a
  mu-coefficient, letting you choose the descent generator.

  If no generator is given, the default descent generator is used.
*/

{
  static CoxWord g(0);

  fprintf(stdout,"first : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  CoxNbr x = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  fprintf(stdout,"second : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (!W->inOrder(x,y)) {
    fprintf(stderr,"the two elements are not in Bruhat order\n");
    return;
  }

  interactive::OutputFile file;
  showMu(file.f(),W->kl(),x,y,W->interface());

  return;
}

void slocus_f ()

/*
  Response to the slocus command. Prints out the singular locus of the
  Schubert variety cl(X_y).
*/

{
  static CoxWord g(0);

  printf("Enter your element (finish with a carriage-return) :\n");

  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputFile file;
  OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),slocusH,traits);
  printSingularLocus(file.f(),y,W->kl(),W->interface(),traits);

  return;
}

void sstratification_f ()

/*
  Response to the slocus command. Prints out the singular locus of the
  Schubert variety cl(X_y).
*/

{
  static CoxWord g(0);

  printf("Enter your element (finish with a carriage-return) :\n");

  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputFile file;
  OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),sstratificationH,traits);
  printSingularStratification(file.f(),y,W->kl(),W->interface(),traits);

  return;
}

void type_f()

/*
  This function sets the type of W, i.e., gets a type and rank from the
  user and sets W to a new group of that type and rank.
*/

{
  CoxGroup* Wloc = interactive::allocCoxGroup();

  if (ERRNO) {
    Error(ERRNO);
  }
  else {
    delete W;
    wgraph_warning = true;
    W = Wloc;
  }

  return;
}

void uneq_f()

/*
  Response to the uneq command.
*/

{
  activate(uneqCommandTree());
  return;
}

namespace uneq {

void klbasis_f()

/*
  Prints out one element in the k-l basis of the group, in the format
  defined by the current output mode.
*/

{
  CoxWord g(0);

  printf("enter your element (finish with a carriage return) :\n");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  uneqkl::HeckeElt h(0);

  W->uneqcBasis(h,y);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  interactive::OutputFile file;
  OutputTraits& traits = W->outputTraits();

  printHeader(file.f(),basisH,traits);
  printAsBasisElt(file.f(),h,W->schubert(),W->interface(),traits);

  return;
}

void lcells_f()

/*
  This function prints out the left cells in the group.
*/

{
  if (!isFiniteType(W)) {
    printFile(stderr,"lcells.mess",MESSAGE_DIR);
    return;
  }

  FiniteCoxGroup* Wf = dynamic_cast<FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillUEMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputFile file;
  OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),lCellsH,traits);
  printLCells(file.f(),Wf->lUneqCell(),Wf->uneqkl(),Wf->interface(),traits);

  return;
}

void lcorder_f()

/*
  Prints the left cell order of the closure of the current group in a file.
  Works only for finite groups.
*/

{
  if (!isFiniteType(W)) {
    printFile(stderr,"lcorder.mess",MESSAGE_DIR);
    return;
  }

  FiniteCoxGroup* Wf = dynamic_cast<FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillUEMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputFile file;
  OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),lCOrderH,traits);
  printLCOrder(file.f(),Wf->uneqkl(),Wf->interface(),traits);

  return;
}

void lrcorder_f()

/*
  Prints the two-sided cell order of the closure of the current group in a
  file. Works only for finite groups.
*/

{
  if (!isFiniteType(W)) {
    printFile(stderr,"uneq/lrcorder.mess",MESSAGE_DIR);
    return;
  }

  FiniteCoxGroup* Wf = dynamic_cast<FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillUEMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputFile file;
  OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),lrCOrderH,traits);
  printLRCOrder(file.f(),Wf->uneqkl(),Wf->interface(),traits);

  return;
}

void lrcells_f()

/*
  This function prints out the two-sided cells in the group. Works only for
  finite groups.
*/

{
  if (!isFiniteType(W)) {
    printFile(stderr,"uneq/lrcells.mess",MESSAGE_DIR);
    return;
  }

  FiniteCoxGroup* Wf = dynamic_cast<FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillUEMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputFile file;
  OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),lrCellsH,traits);
  printLRCells(file.f(),Wf->lrUneqCell(),Wf->uneqkl(),Wf->interface(),traits);

  return;
}

void mu_f()

/*
  Response to the "mu" command. Prints out a single mu-coefficient. When
  the generator s act on the left, we use the fact that mu(left_s,x,y) is
  equal to mu(right_s,x^{-1},y^{-1}) to go over to the right action. This
  is necessary because the mu-tables are kept only for right actions.
*/

{
  static CoxWord g(0);
  bool leftAction = false;

  fprintf(stdout,"generator : ");
  Generator s = getGenerator(W);

  if (s >= W->rank()) { // action is on the left
    s -= W->rank();
    leftAction = true;
  }

  fprintf(stdout,"first : ");
  g = interactive::getCoxWord(W);
  if (leftAction)
    W->inverse(g);
  if (!W->isDescent(g,s)) { // mu(s,x,y) is undefined
    fprintf(stderr,"xs is greater than x\n");
    return;
  }
  CoxNbr x = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  fprintf(stdout,"second : ");
  g = interactive::getCoxWord(W);
  if (leftAction)
    W->inverse(g);
  if (W->isDescent(g,s)) { // mu(s,x,y) is undefined
    fprintf(stderr,"ys is smaller than y\n");
    return;
  }
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (x == y) {
    fprintf(stderr,"the two elements are equal\n");
    return;
  }

  if (!W->inOrder(x,y)) {
    fprintf(stderr,"the two elements are not in Bruhat order\n");
    return;
  }

  const uneqkl::MuPol& mu = W->uneqmu(s,x,y);
  if (ERRNO) {
    Error(ERRNO,x,y);
    return;
  }

  print(stdout,mu,"v");
  printf("\n");

  return;
}

void pol_f()

/*
  Response to the "pol" command. Prints out a single polynomial.
*/

{
  static CoxWord g(0);

  fprintf(stdout,"first : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  CoxNbr x = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  fprintf(stdout,"second : ");
  g = interactive::getCoxWord(W);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }
  CoxNbr y = W->extendContext(g);
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  if (!W->inOrder(x,y)) {
    fprintf(stderr,"the two elements are not in Bruhat order\n");
    return;
  }

  const uneqkl::KLPol& pol = W->uneqklPol(x,y);
  if (ERRNO) {
    Error(ERRNO,x,y);
    return;
  }

  print(stdout,pol,"q");
  printf("\n");

  return;
}

void rcells_f()

/*
  This function prints out the left cells in the group. Works only for finite
  groups.
*/

{
  if (!isFiniteType(W)) {
    printFile(stderr,"rcells.mess",MESSAGE_DIR);
    return;
  }

  FiniteCoxGroup* Wf = dynamic_cast<FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillUEMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputFile file;
  OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),rCellsH,traits);
  printRCells(file.f(),Wf->rUneqCell(),Wf->uneqkl(),Wf->interface(),traits);

  return;
}

void rcorder_f()

/*
  Prints the right cell order of the current group in a file. Works only for
  finite groups.
*/

{
  if (!isFiniteType(W)) {
    printFile(stderr,"rcorder.mess",MESSAGE_DIR);
    return;
  }

  FiniteCoxGroup* Wf = dynamic_cast<FiniteCoxGroup*> (W);

  Wf->fullContext();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  Wf->fillUEMu();
  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  OutputFile file;
  OutputTraits& traits = Wf->outputTraits();

  printHeader(file.f(),rCOrderH,traits);
  printRCOrder(file.f(),Wf->uneqkl(),Wf->interface(),traits);

  return;
}

};

};

namespace commands {

void interface::abort_f()

/*
  Aborts the interface modification. Bypasses the exit function,
  so that there is no further checking of the abandoned choices.
*/

{
  delete in_buf;
  in_buf = 0;
  treeStack.pop();

  return;
}

void interface::alphabetic_f()

/*
  Sets i/o to the standard alphabetic conventions : the symbols are the
  alphabetic sequence, prefix and postfix are empty, and the separator is
  "." iff the rank is > 26.
*/

{
  delete in_buf;
  in_buf = new GroupEltInterface(W->rank(),Alphabetic());
  W->interface().setIn(*in_buf);
  W->interface().setOut(*in_buf);

  return;
}

void interface::bourbaki_f()

/*
  Sets Bourbaki conventions. This means that the ordering is reverted in
  types B and D, and symbols as well.

  NOTE : currently not implemented for affine groups.
*/

{
  delete in_buf;
  in_buf = new GroupEltInterface(W->interface().inInterface());
  in::bourbaki_f();
  W->interface().setIn(*in_buf);

  delete in_buf;
  in_buf = new GroupEltInterface(W->interface().outInterface());
  out::bourbaki_f();
  W->interface().setOut(*in_buf);

  return;
}

void interface::decimal_f()

/*
  Sets i/o to the standard decimal conventions : the symbols are the
  decimal sequence, prefix and postfix are empty, and the separator is
  "." iff the rank is > 9.
*/

{
  delete in_buf;
  in_buf = new GroupEltInterface(W->rank(),Decimal());
  W->interface().setIn(*in_buf);
  W->interface().setOut(*in_buf);

  return;
}

void interface::default_f()

/*
  Sets i/o settings to the default style. This means that we use decimal
  symbols, no prefix or postfix, and separator "." only when rank is >= 10.
  The ordering is the internal default ordering.
*/

{
  delete in_buf;
  in_buf = new GroupEltInterface(W->rank());

  W->interface().setIn(*in_buf);
  W->interface().setOut(*in_buf);

  W->interface().setOrder(identityOrder(W->rank()));
  W->interface().setDescent(Default());
  W->setOutputTraits(Pretty());

  return;
}

void interface::gap_f()

/*
  Sets i/o settings to GAP style. This means first of all that Bourbaki
  conventions are adopted; decimal symbols are used for i/o with prefix
  "[", separator "," and postfix "]". Furthermore, output to files is
  done in GAP style, which produces files that are directly legible by
  GAP3.
*/

{
  delete in_buf;
  in_buf = new GroupEltInterface(W->rank(),GAP());

  in::bourbaki_f();
  W->interface().setIn(*in_buf);

  out::bourbaki_f();
  W->interface().setOut(*in_buf);
  W->interface().setDescent(GAP());
  W->setOutputTraits(GAP());

  return;
}

void interface::hexadecimal_f()

/*
  Sets i/o to the standard hexadecimal conventions : the symbols are the
  hexadecimal sequence, prefix and postfix are empty, and the separator is
  "." iff the rank is > 15.
*/

{
  delete in_buf;
  in_buf = new GroupEltInterface(W->rank(),Hexadecimal());
  W->interface().setIn(*in_buf);

  W->interface().setOut(*in_buf);

  return;
}

void interface::in_f()

{
  activate(inCommandTree());
  return;
}

void interface::ordering_f()

/*
  Allows the user to change the generator ordering.
*/

{
  static Permutation in_order(W->rank());

  changeOrdering(W,in_order);

  if (ERRNO) {
    Error(ERRNO);
    return;
  }

  W->setOrdering(in_order);

  return;
}

void interface::out_f()

{
  activate(outCommandTree());
  return;
}

void interface::permutation_f()

/*
  Activates permutation i/o.
*/

{
  using namespace coxeter;

  if (!isTypeA(W->type())) {
    printFile(stderr,"permutation.mess",MESSAGE_DIR);
    return;
  }

  TypeACoxGroup* WA = dynamic_cast<TypeACoxGroup*>(W);

  WA->setPermutationInput(true);
  WA->setPermutationOutput(true);

  W->interface().setOrder(identityOrder(W->rank()));
  W->interface().setDescent(Default());
  W->setOutputTraits(Pretty());

  return;
}

/**
  Resets a symbol in in_buf (this will become either an input or an output
  symbol).
*/
void interface::symbol_f()
{
  static String buf(0);

  const Interface& I = W->interface();
  Generator s = undef_generator;
  reset(buf);

  do {
    if (ERRNO)
      Error(ERRNO);
    printf("enter the generator symbol you wish to change, ? to abort:\n");
    getInput(stdin,buf,0);
    if (buf[0] == '?')
      return;
    io::skipSpaces(buf,0);
    Token tok;
    I.symbolTree().find(buf,0,tok);
    if (tokenType(tok) != generator_type)/* error */
      ERRNO = NOT_GENERATOR;
    else
      s = tok-1;
  } while (ERRNO);

  printf("enter the new symbol (finish with a carriage return):\n");
  getInput(stdin,buf,0);
  in_buf->setSymbol(s,buf);

  return;
}

void interface::terse_f()

/*
  Sets i/o settings to terse style. This style is meant for outputting files
  that are easily parsed by computer. Th ordering of the generators is left
  untouched, and can be set independently by the user. Decimal output symbols
  are chosen, prefix is set to "[", postfix to "]" and separator to ",".
*/

{
  delete in_buf;
  in_buf = new GroupEltInterface(W->rank(),GAP());
  W->interface().setIn(*in_buf);
  W->interface().setOut(*in_buf);

  W->interface().setDescent(Default());
  W->setOutputTraits(Terse());

  return;
}

void interface::in::alphabetic_f()

/*
  Sets the input symbols to the alphabetic sequence.
*/

{
  const String* alpha = alphabeticSymbols(in_buf->symbol.size());

  for (Ulong j = 0; j < in_buf->symbol.size(); ++j) {
    in_buf->symbol[j] = alpha[j];
  }

  return;
}

void interface::in::bourbaki_f()

/*
  Sets Bourbaki conventions for input. This means reverting the ordering
  of the input symbols in types B and D.
*/

{
  const Type& x = W->type();

  if (!isFiniteType(x))
    return;
  if (!(isTypeB(x) || isTypeD(x)))
    return;

  for (Generator s = 0; s < W->rank(); ++s) {
    in_buf->symbol[s] = W->interface().inSymbol(W->rank()-s-1);
  }

  return;
}

void interface::in::decimal_f()

/*
  Sets the input symbols to the decimal sequence.
*/

{
  const String* dec = decimalSymbols(in_buf->symbol.size());

  for (Ulong j = 0; j < in_buf->symbol.size(); ++j) {
    in_buf->symbol[j] = dec[j];
  }

  return;
}

void interface::in::default_f()

/*
  Sets the input interface to the default style.
*/

{
  delete in_buf;
  in_buf = new GroupEltInterface(W->rank());

  return;
}

void interface::in::gap_f()

/*
  Sets the input interface to GAP style, and enforces Bourbaki conventions.
*/

{
  delete in_buf;
  in_buf = new GroupEltInterface(W->rank(),GAP());
  in::bourbaki_f();

  return;
}

void interface::in::hexadecimal_f()

/*
  Sets the input symbols to the hexadecimal sequence.
*/

{
  const String* hex = hexSymbols(in_buf->symbol.size());

  for (Ulong j = 0; j < in_buf->symbol.size(); ++j) {
    in_buf->symbol[j] = hex[j];
  }

  return;
}

void interface::in::permutation_f()

/*
  Sets input to permutation mode (type A only.)
*/

{
  using namespace coxeter;

  if (!isTypeA(W->type())) {
    printFile(stderr,"permutation.mess",MESSAGE_DIR);
    return;
  }

  TypeACoxGroup* WA = dynamic_cast<TypeACoxGroup*>(W);
  WA->setPermutationInput(true);

  delete in_buf;
  in_buf = 0;

  return;
}

void interface::in::postfix_f()

/*
  Resets the input postfix.
*/

{
  printf("Enter the new input postfix (finish with a carriage return):\n");
  String buf(0);
  getInput(stdin,buf,0);
  in_buf->setPostfix(buf);
  return;
}

void interface::in::prefix_f()

/*
  Resets the input prefix.
*/

{
  printf("Enter the new input prefix (finish with a carriage return):\n");
  String buf(0);
  getInput(stdin,buf,0);
  in_buf->setPrefix(buf);
  return;
}

void interface::in::separator_f()

/*
  Resets the input separator.
*/

{
  printf("Enter the new input separator (finish with a carriage return):\n");
  String buf(0);
  getInput(stdin,buf,0);
  in_buf->setSeparator(buf);
  return;
}

void interface::in::terse_f()

/*
  Sets the input interface to terse style (the same as GAP style).
*/

{
  delete in_buf;
  in_buf = new GroupEltInterface(W->rank(),GAP());

  return;
}

void interface::out::alphabetic_f()

/*
  Changes the output symbols to hexadecimal.
*/

{
  const String* alpha = alphabeticSymbols(in_buf->symbol.size());

  for (Ulong j = 0; j < in_buf->symbol.size(); ++j) {
    in_buf->symbol[j] = alpha[j];
  }

  return;
}

void interface::out::bourbaki_f()

/*
  Sets Bourbaki conventions for input. This means reverting the ordering
  of the output symbols in types B and D, and setting the output ordering
  to the reverse of the default one.
*/

{
  const Type& x = W->type();

  if (!isFiniteType(x))
    return;
  if (!(isTypeB(x) || isTypeD(x))) {
    W->setOrdering(identityOrder(W->rank()));
    return;
  }

  for (Generator s = 0; s < W->rank(); ++s) {
    in_buf->symbol[s] = W->interface().outSymbol(W->rank()-s-1);
  }

  Permutation a(W->rank());

  for (Generator s = 0; s < W->rank(); ++s) {
    a[s] = W->rank()-1-s;
  }

  W->setOrdering(a);

  return;
}

void interface::out::decimal_f()

/*
  Changes the output symbols to decimal.
*/

{
  const String* dec = decimalSymbols(in_buf->symbol.size());

  for (Ulong j = 0; j < in_buf->symbol.size(); ++j) {
    in_buf->symbol[j] = dec[j];
  }

  return;
}

void interface::out::default_f()

/*
  Sets output styles to the default style.
*/

{
  delete in_buf;
  in_buf = new GroupEltInterface(W->rank());
  W->setOrdering(identityOrder(W->rank()));

  W->setOutputTraits(Pretty());

  return;
}

void interface::out::gap_f()

/*
  Sets output styles to GAP style, and enforces Bourbaki conventions.
*/

{
  delete in_buf;
  in_buf = new GroupEltInterface(W->rank(),GAP());
  W->setOrdering(identityOrder(W->rank()));
  out::bourbaki_f();

  W->interface().setDescent(GAP());
  W->interface().setOut(*in_buf); // has to be done here so that output traits
                                  // will be correct.
  W->setOutputTraits(GAP());

  return;
}

void interface::out::hexadecimal_f()

/*
  Changes the output symbols to hexadecimal.
*/

{
  const String* hex = hexSymbols(in_buf->symbol.size());

  for (Ulong j = 0; j < in_buf->symbol.size(); ++j) {
    in_buf->symbol[j] = hex[j];
  }

  return;
}

void interface::out::permutation_f()

/*
  Sets output to permutation mode (type A only.)
*/

{
  using namespace coxeter;

  if (!isTypeA(W->type())) {
    printFile(stderr,"permutation.mess",MESSAGE_DIR);
    return;
  }

  TypeACoxGroup* WA = dynamic_cast<TypeACoxGroup*>(W);

  WA->setPermutationOutput(true);

  W->interface().setOrder(identityOrder(W->rank()));
  W->interface().setDescent(Default());
  W->setOutputTraits(Pretty());

  delete in_buf;
  in_buf = 0;

  return;
}

void interface::out::postfix_f()

/*
  Resets the output postfix.
*/

{
  printf("enter the new output postfix (finish with a carriage return):\n");
  String buf(0);
  getInput(stdin,buf,0);
  in_buf->setPostfix(buf);
  return;
}

void interface::out::prefix_f()

/*
  Resets the output prefix.
*/

{
  printf("Enter the new output prefix (finish with a carriage return):\n");
  String buf(0);
  getInput(stdin,buf,0);
  in_buf->setPrefix(buf);
  return;
}

void interface::out::separator_f()

/*
  Resets the output separator.
*/

{
  printf("Enter the new output separator (finish with a carriage return):\n");
  String buf(0);
  getInput(stdin,buf,0);
  in_buf->setSeparator(buf);
  return;
}

void interface::out::terse_f()

/*
  Sets output styles to terse style (the same as GAP style).
*/

{
  delete in_buf;
  in_buf = new GroupEltInterface(W->rank(),GAP());

  W->interface().setDescent(Default());
  W->interface().setOut(*in_buf); // has to be done here so that output
                                  // traits will be correct
  W->setOutputTraits(Terse());

  return;
}

};

/*****************************************************************************

        Chapter VI -- Miscellaneous.

  This section contains some auxiliary functions :

  - printCommands(file,tree) : prints info about the various commands
    on the tree;
  - printCommandTree(file,cell) : the recursive function doing the actual
    work;
  - interface_entry() : entry function for the interface mode;
  - interface_exit() : exit function for the interface mode;
  - interface::in_entry() : entry function for the interface::in mode;
  - interface::in_exit() : exit function for the interface::in mode;
  - interface::out_entry() : entry function for the interface::out mode;
  - interface::out_exit() : exit function for the interface::out mode;
  - main_entry() : entry function for the main mode;
  - main_exit() : exit function for the main mode;
  - uneq_entry() : entry function for the uneq mode;
  - uneq_exit() : exit function for the uneq mode;

 *****************************************************************************/

namespace commands {

void printCommands(FILE* file, CommandTree* tree)

/*
  Prints one line for each command on the tree (sorted in alphabetical order)
  with the name of the command and the information contained in the tag field.
*/

{
  printCommandTree(file,tree->root()->left);
  return;
}

};

CoxGroup* commands::currentGroup()

/*
  Returns the "current" Coxeter group.

  NOTE : this will probably have to be refined in the future.
*/

{
  return W;
}

namespace {

void printCommandTree(FILE* file, DictCell<CommandData>* cell)

{
  if (cell == 0)
    return;

  if (cell->fullname) { /* print command info */
    CommandData* cd = cell->value();
    fprintf(file,"  - %s : %s;\n",cd->name.ptr(),cd->tag.ptr());
  };

  printCommandTree(file,cell->left);
  printCommandTree(file,cell->right);

  return;
}

void interface_entry()

{
  commands::interface::in_buf = new GroupEltInterface(W->rank());
  return;
}

void interface_exit()

{
  delete commands::interface::in_buf;
  commands::interface::in_buf = 0;
  return;
}

void main_entry()

/*
  Sets the type and rank. This is used as entry function for the main mode,
  and as the function that restarts the program when we change the type.

  NOTE : error handling should be done by the calling function.

  NOTE : something should be done about deallocating before reallocating!
*/

{
  W = interactive::allocCoxGroup();

  /* an error may be set here */

  return;
}

};

namespace commands {

void interface::in_entry()

/*
  Entry function to the input interface modification mode. The global variable
  in_buf is originally set to value for the current group.
*/

{
  Permutation a(W->interface().order());
  a.inverse();

  printf("current input symbols are the following :\n\n");
  printInterface(stdout,W->interface().inInterface(),a);
  printf("\n");

  in_buf = new GroupEltInterface(W->interface().inInterface());

  return;
}

void interface::in_exit()

/*
  Exit function from the input modification mode. It checks if the
  modifications made by the user are consistent with the ones that will
  remain from the old interface; if yes, it confirms the modifications
  and exits peacefully; if not, it prints out the problems and keeps
  the user in the mode.
*/

{
  if (in_buf == 0) // hack to prevent execution in special cases
    return;

  Permutation a(W->interface().order());
  a.inverse();

  /* at this point in_buf holds the full putative new interface; we
   need to check for reserved or repeated non-empty symbols */

  const String* str = checkLeadingWhite(*in_buf);

  if (str) {
    Error(LEADING_WHITESPACE,in_buf,&W->interface().inInterface(),&a,str);
    goto error_exit;
  }

  str = checkReserved(*in_buf,W->interface());

  if (str) {
    Error(RESERVED_SYMBOL,in_buf,&W->interface().inInterface(),&a,str);
    goto error_exit;
  }

  if (!checkRepeated(*in_buf)) {
    Error(REPEATED_SYMBOL,in_buf,&W->interface().inInterface(),&a);
    goto error_exit;
  }

  /* if we reach this point, the new interface is ok */

  printf("new input symbols:\n\n");
  printInterface(stdout,*in_buf,a);
  printf("\n");

  W->interface().setIn(*in_buf);

  return;

 error_exit:
  ERRNO = ERROR_WARNING;
  return;
}

void interface::out_entry()

/*
  Entry function to the output interface modification mode. The global variable
  in_buf is originally set to value for the current group.
*/

{
  delete in_buf;
  in_buf = new GroupEltInterface(W->interface().outInterface());

  Permutation a(W->interface().order());
  a.inverse();

  printf("current output symbols are the following :\n\n");
  printInterface(stdout,*in_buf,W->interface().inInterface(),a);
  printf("\n");

  return;
}

void interface::out_exit()

/*
  Exit function for the output modification mode. No checking is necessary
  here.
*/

{
  if (in_buf == 0) // hack to prevent execution in special cases
    return;

  Permutation a(W->interface().order());
  a.inverse();

  printf("new output symbols:\n\n");
  printInterface(stdout,*in_buf,W->interface().inInterface(),a);
  printf("\n");

  W->interface().setOut(*in_buf);

  return;
}

};

namespace {

void main_exit()

/*
  Symmetric function to main_entry. Should undo what main_entry did.
*/

{
  delete W;
  wgraph_warning = true;
  return;
}

void uneq_entry()

{
  W->activateUEKL();
  return;
}

void uneq_exit()

/*
  We keep the unequal-parameter context, because one might want to go back
  and forth between the unequal and the ordinary case.
*/

{
  return;
}

};
