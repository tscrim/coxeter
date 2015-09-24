/*
  This is commands.h
  
  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file main.cpp for full copyright notice
*/

#ifndef COMMANDS_H  /* guard against multiple inclusions */
#define COMMANDS_H

#include "globals.h"

namespace commands {
  using namespace globals;
};

#include "coxgroup.h"

/******** type declarations ************************************************/

namespace commands {
  struct CommandData;
  class CommandTree;
};

/******** constants ********************************************************/

namespace commands {
  extern void (*default_help)();
};

/******** function declarations ********************************************/

namespace commands {
  coxgroup::CoxGroup* currentGroup();
  void default_error(char* str);
  void execute();
  CommandTree* interfaceCommandTree();
  CommandTree* mainCommandTree();
  void printCommands(FILE* file, CommandTree* tree);
  void relax_f();
  void run();
  CommandTree* uneqCommandTree();
  namespace interface {
    CommandTree* inCommandTree();
    CommandTree* outCommandTree();
  };
};

/******** Type definitions *************************************************/

#include "dictionary.h"
#include "io.h"

namespace commands {
  using namespace dictionary;
  using namespace io;
};

namespace commands {

struct CommandData {
  String name;
  String tag;
  void (*action)();
  void (*help)();
  bool autorepeat;
/* Constructors and destructors */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(CommandData));}
  CommandData(const char* const& str, const char* const& t, void (*a)(), 
	      void (*h)(), bool rep);
  ~CommandData();
};

class CommandTree:public Dictionary<CommandData> {
 private:
  String d_prompt;
  CommandTree* d_help;
  void (*d_entry)();
  void (*d_error)(char* str);
  void (*d_exit)();
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(CommandTree));}
  CommandTree(const char *str, void (*action)(), void (*entry)() = &relax_f,
	      void (*error)(char*) = &default_error, 
	      void (*exit)() = &relax_f, void (*h)() = 0);
  ~CommandTree();
/* modifiers */
  void add(const char* name, const char* tag, void (*action)(), 
	   void (*help)() = default_help, bool rep = true);
  void setAction(const char* str, void (*a)());
  void setRepeat(const char* str, bool b);
  void setEntry(void (*a)());                                    /* inlined */
/* accessors */
  void prompt() const;
  void entry() const;                                            /* inlined */
  void error(char *str) const;                                   /* inlined */
  void exit() const;                                             /* inlined */
  CommandTree* helpMode() const;                                 /* inlined */
};

};

/******** Inline definitions *********************************************/

namespace commands {

inline void CommandTree::setEntry(void (*a)()) {d_entry = a;}
inline void CommandTree::entry() const {return d_entry();}
inline void CommandTree::error(char *str) const {return d_error(str);}
inline void CommandTree::exit() const {return d_exit();}
inline CommandTree* CommandTree::helpMode() const {return d_help;}

};

#endif
