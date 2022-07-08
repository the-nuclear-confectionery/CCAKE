#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

#include "../include/formatted_output.h"

using std::ceil;
using std::cout;
using std::endl;
using std::floor;
using std::istringstream;
using std::string;
using std::to_string;
using std::vector;

namespace formatted_output
{
  //----------------------------------------------------------------------------
  void print_bullet( string message, string b, int indent )
  {
    int max_length = 80;
    int bn = b.length();

    string bpad = string(bn, ' ');
    string indentstring = string(indent, ' ');

    if ( message.length() < max_length - indent - bn )
      cout << indentstring << b << message;
    else
    {
      istringstream iss(message);
      string word = "";
      int current_length = indent + bn;
      cout << indentstring << b;
      while (iss >> word)
      {
        if (current_length + word.length() > max_length)
        {
          cout << "\n" << indentstring << bpad << " " << word << " ";
          current_length = word.length() + indent + bn + 1;
        }
        else
        {
          cout << word << " " ;
          current_length += word.length() + 1;
        }
      }
    }

    cout << endl;

    return;
  }

  //----------------------------------------------------------------------------
  void print_flexible( string message, char lc, char rc )
  {
    bool RHbuffer = true;
    int max_length = 80;
    cout << string(max_length/2, lc) + string(max_length/2, rc) << "\n";
    int buffer = ( RHbuffer ) ? 2 : 0;
    max_length -= buffer;
    if ( message.length() <= max_length - 4 )
      cout << lc << " " << message << " " << rc << "\n";
    else
    {
      istringstream iss(message);
      string word = "";
      int current_length = 2;
      cout << lc << " ";
      while (iss >> word)
      {
        if (current_length + word.length() > max_length)
        {
          if (RHbuffer)
            cout << std::string( max_length - current_length + 1, ' ' ) << rc;
          cout << "\n" << lc << " " << word << " ";
          current_length = word.length() + 3;
        }
        else
        {
          cout << word << " " ;
          current_length += word.length() + 1;
        }
      }
      if (RHbuffer)
        cout << std::string( max_length - current_length + 1, ' ' ) << rc;
    }

    max_length += buffer;
    cout << "\n" << string(max_length/2, lc) + string(max_length/2, rc) << endl;

    return;
  }

  //============================================================================
  //----------------------------------------------------------------------------
  void announce( string message )
  {
    cout << "\n\n";
    int max_length = 80;
    string wrapper = string(max_length/2, '<') + string(max_length/2, '>');
    string padded_message = string(6, ' ') + message + string(6, ' ');

    if ( padded_message.length() <= max_length - 2 )
    {
      cout << wrapper << "\n";
      double half_width = 0.5*(max_length - padded_message.length());
      string left_padding  = string( (int)floor(half_width), '<' );
      string right_padding = string( (int)ceil(half_width), '>' );
      cout << left_padding << padded_message << right_padding << "\n";
      cout << wrapper << "\n";
    }
    else
      print_flexible( message, '<', '>' );

    cout << endl;

    return;
  }

  //----------------------------------------------------------------------------
  void summarize(string message)
  {
    cout << endl;
    print_bullet( message, "", 0 );
  }

  //----------------------------------------------------------------------------
  void report(string message)
  {
    cout << endl;
    print_bullet( message, "- ", 2 );
  }

  //----------------------------------------------------------------------------
  void update(string message)
  {
    cout << endl;
    print_bullet( message, "* ", 4 );
  }


  //----------------------------------------------------------------------------
  void detail(string message)
  {
    print_bullet( message, "--> ", 6 );
  }


  //----------------------------------------------------------------------------
  void comment(string message)
  {
    print_bullet( message, "+ ", 10 );
  }

  
}