/**********************************************************************************************
    
    this little program read events from a lhe file and randomnly send intermediate/final 
    particles in their antiparticles.
    this operation remove a possible charge asymmetry produced by a bug in madgraph
     
    compile with ---> c++ -std=c++0x -o flavourFlip flavourFlip.cpp

**********************************************************************************************/

#include "LHEF.h"
#include <iomanip>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <cmath>
#include <random>
#include <functional>

using namespace std ;

//*********************************************************************************************

int main(int argc, char ** argv) 
{
    if(argc < 3)
    {
        cout << "Usage:   " << argv[0] 
             << " input.lhe output.lhe" << endl ;
        return -1;
    }

    std::ifstream ifs (argv[1]) ;
    LHEF::Reader reader (ifs) ;

    ofstream outputStream (argv[2]) ;
    LHEF::Writer writer (outputStream) ;

    writer.headerBlock() << reader.headerBlock ;
    writer.initComments() << reader.initComments ;
    writer.heprup = reader.heprup ;
    writer.init () ;

    default_random_engine generator;
    uniform_int_distribution<int> distribution(0,1);
    auto coin = bind (distribution, generator);

    int iEvent=0;
    while (reader.readEvent ()) 
    {
        iEvent++;
        if ( iEvent % 10000 == 0 )
        {
            cout << iEvent << "  processed" << endl;
        }
        if ( reader.outsideBlock.length() ) 
        {
            std::cout << reader.outsideBlock;
        }
        int coin_tmp = coin();
        if (coin_tmp == 1)
        {
            // loop over particles in the event
            for (int iPart = 2 ; iPart < reader.hepeup.IDUP.size (); iPart++) 
            {
                // particle status (incoming==-1, intermediate==2, outgoing==1)
                int pSTATUS = reader.hepeup.ISTUP.at (iPart);
                int pPID = reader.hepeup.IDUP.at (iPart);          
                int mother1 = (reader.hepeup.MOTHUP.at (iPart)).first;
                int mother1PID = reader.hepeup.IDUP.at (mother1-1);   
                int mother2 = (reader.hepeup.MOTHUP.at (iPart)).second;
                int mother2PID = reader.hepeup.IDUP.at (mother2-1);   
                if ( (pSTATUS == 1 && abs(mother1PID) == 24) || (pSTATUS == 2 && abs(pPID) == 24) )
                {
                    reader.hepeup.IDUP.at (iPart) = -reader.hepeup.IDUP.at (iPart);
                    if( pPID < 7 )
                    {
                        int colour = (reader.hepeup.ICOLUP.at (iPart)).first;
                        int anticolour = (reader.hepeup.ICOLUP.at (iPart)).second;
                        (reader.hepeup.ICOLUP.at (iPart)).first = anticolor;
                        (reader.hepeup.ICOLUP.at (iPart)).second = color;
                    }         
                }
            }
        } 
        writer.eventComments() << reader.eventComments;
        writer.hepeup = reader.hepeup;
        writer.writeEvent();
    } 
    
    return 0 ;
}
