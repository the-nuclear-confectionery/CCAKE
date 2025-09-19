#ifndef WELCOME_H
#define WELCOME_H

#include <string>

#include "config.h"

using std::cout;
using std::endl;

namespace message
{

void welcome()
  {
    cout << "================================================================================\n"
         << "=                          ____ ____    _    _  _______                        =\n"
            "=                         / ___/ ___|  / \\  | |/ / ____|                       =\n"
            "=                        | |  | |     / _ \\ | ' /|  _|                         =\n"
            "=                        | |__| |___ / ___ \\| . \\| |___                        =\n"
            "=                         \\____\\____/_/   \\_\\_|\\_\\_____|                       =\n"
            "=                                                                              =\n"
         << "================================================================================\n"
         << "= Version: " << ccake_VERSION_MAJOR << "." << ccake_VERSION_MINOR << "\n"
         << "================================================================================\n";
    cout << "= Authors:\n";
    cout << "= Christopher Plumberg\n";
    cout << "= Travis Dore\n";
    cout << "= Dekrayat Almaalol\n";
    cout << "= Débora Mroczek\n";
    cout << "= Jaquelyn Noronha-Hostler\n";
    cout << "= Lydia Spychalla\n";
    cout << "= Jordi Salinas San Martin\n";
    cout << "= Willian M. Serenone\n";
    cout << "= Kevin P. Pala\n";
    cout << "= Surkhab Kaur Virk\n";
    cout << "= Isabella Danhoni\n";
    cout << "= Isaac Long\n";
    cout << "================================================================================" << endl;

  }

  void guilty_welcome()
  {
    cout << endl << endl;
    cout << "-------------------------------------------------------------------------------------"
         << endl
         << "______ _     _                            _ _                                _  ___  " << endl
         << "|  _  (_)   | |                          (_) |                              | ||__ \\ " << endl
         << "| | | |_  __| |  _   _  ___  _   _    ___ _| |_ ___   _   _ ___   _   _  ___| |_  ) |" << endl
         << "| | | | |/ _` | | | | |/ _ \\| | | |  / __| | __/ _ \\ | | | / __| | | | |/ _ \\ __|/ / " << endl
         << "| |/ /| | (_| | | |_| | (_) | |_| | | (__| | ||  __/ | |_| \\__ \\ | |_| |  __/ |_|_|  " << endl
         << "|___/ |_|\\__,_|  \\__, |\\___/ \\__,_|  \\___|_|\\__\\___|  \\__,_|___/  \\__, |\\___|\\__(_)  " << endl
         << "                  __/ |                                            __/ |             " << endl
         << "                 |___/                                            |___/              " << endl
         << endl;
    cout << "-------------------------------------------------------------------------------------"
         << endl;
    cout << "Authors: " << endl;
    cout << "Christopher Plumberg" << endl;
    cout << "Travis Dore" << endl;
    cout << "Dekrayat Almaalol" << endl;
    cout << "Jaquelyn Noronha-Hostler" << endl;
    cout << "Lydia Spychalla" << endl;
    cout << "Jordi Salinas San Martin" << endl ;
    cout << "Willian M. Serenone"<< endl << endl ;
    cout << "Kevin P. Pala" << endl;
    cout << "Surkhab Kaur Virk" << endl;
    cout << "Isabella Danhoni" << endl;
    cout << "Isaac Long" << endl;
    cout << "-------------------------------------------------------------------------------------"
         << endl << endl;

  }

}

#endif