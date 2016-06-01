// This file is part of Swept2D
// Copyright (C) 2015 Qiqi Wang, qiqi@mit.edu AND Maitham Alhubail, hubailmm@mit.edu
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU Affero General Public License as
// published by the Free Software Foundation, either version 3 of the
// License, or (at your option) an later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT An WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU Affero General Public License for more details.
//
// You should have received a copy of the GNU Affero General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef H_MPIGLOBALS2D
#define H_MPIGLOBALS2D

#include <mpi.h>
#include <algorithm>
#include <string>
#include "ProcessGraph.h"
#include "ServerSocket.h"
#include "SocketException.h"
#include <string>
#include <iostream>
#include <vector>
#include <cctype>
#include "base64.h"
#ifdef _WIN32
#include <Windows.h>
#else
#include <pthread.h>
#endif

ProcessGraph pg;
string constantsArray;
string resultsArray;
vector<string> remoteCommands;
vector<int> remoteCodes;
Socket *connection;
unsigned char* constantsArrayBytes;

#ifdef _WIN32
HANDLE threadHd;
DWORD  threadId;
HANDLE remoteCommandsThreadMutexLock;
#else
pthread_t remoteCommandsThread;
pthread_mutex_t remoteCommandsThreadMutexLock;
#endif

static inline std::string &ltrim(std::string &s) 
{
	s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

static inline std::string &rtrim(std::string &s) 
{
	s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
    return s;
}

static inline std::string &trim(std::string &s) 
{
	return ltrim(rtrim(s));
}

int ijToIndex(int i,int j,int xNodes,int yNodes)
{
	int ii,jj;
	ii = i;
	jj = j;
	if(i==-1)      ii = xNodes-1;
	else if(i==xNodes) ii = 0;
	if(j==-1)      jj = yNodes-1;
	else if(j==yNodes) jj = 0;

	int index = jj * xNodes + ii;
	return index;
}

ProcessGraph getProcessGraph(int myrank,int xNodes,int yNodes,int size)
{
	ProcessGraph pg;
	pg.rank = myrank;
	int j = (myrank % (xNodes*yNodes)) / xNodes;
	int i = myrank % xNodes;
	
	pg.iIndex  = i;
	pg.jIndex  = j;
	pg.xNodes  = xNodes;
	pg.yNodes  = yNodes;
	pg.mpiSize = size;
	pg.Wrank   = ijToIndex(i-1,j,xNodes,yNodes);
	pg.Erank   = ijToIndex(i+1,j,xNodes,yNodes);
	pg.Srank   = ijToIndex(i,j+1,xNodes,yNodes);
	pg.Nrank   = ijToIndex(i,j-1,xNodes,yNodes);				
	return pg;
}

void sendCommandToSolver(int code,string &data)
{
	#ifdef _WIN32
	WaitForSingleObject(remoteCommandsThreadMutexLock,INFINITE);
	#else
	pthread_mutex_lock(&remoteCommandsThreadMutexLock);
	#endif
	remoteCommands.push_back(data);
	remoteCodes.push_back(code);
	#ifdef _WIN32
	ReleaseMutex(remoteCommandsThreadMutexLock);
	#else
	pthread_mutex_unlock(&remoteCommandsThreadMutexLock);
	#endif
	code = 0;
	data.clear();
}

void *listenToCommands(void *ptr)
{
	ServerSocket server(30000);
	string welcome("Connected to Swept2D Solver.");
	char messageLength[9];

	while(true)
	{
		try
		{
			string message;
			connection = server.Accept();
			printf("Client Connected!!\n");		
			sprintf(messageLength,"%8u",(int)welcome.size());
			message.append(messageLength);
			message.append(welcome);
			connection->sendLine(message);
			while(true)
			{
				string data;
				string response;
				int code;
				if(connection->recvLine(response) == 0)
				{
					std::replace( response.begin(), response.end(), '\r', ' ');
					std::replace( response.begin(), response.end(), '\n', ' ');
					trim(response);
					std::transform(response.begin(), response.end(),response.begin(), ::toupper);
					//printf("Got: %s\n",response.c_str());
					if(response.compare("UPDATE CONSTANTS") == 0)
					{
						string message;
						response.clear();
						response = "OK! Please send the new values for the constants as BYTES Base64 encoded!";
						sprintf(messageLength,"%8u",(int)response.size());
						message.append(messageLength);
						message.append(response);
						//connection->sendLine(message);
						constantsArray.clear();
						connection->recvLine(constantsArray);
						std::replace( constantsArray.begin(), constantsArray.end(), '\r', ' ');
						std::replace( constantsArray.begin(), constantsArray.end(), '\n', ' ');
						trim(constantsArray);
						response.clear();
						response += "Data Received! Will apply the changes ASAP!";
						sprintf(messageLength,"%8u",(int)response.size());
						message.clear();
						message.append(messageLength);
						message.append(response);
						//connection->sendLine(message);
						int outSize = 1024000000;
						base64decode((char*)constantsArray.c_str(),constantsArray.size(),constantsArrayBytes,(size_t*)&outSize);
						printf("Received Constants Update Data.  Total Bytes: %d\n",outSize);
						data = "UPDATE";
						code = 1;
						sendCommandToSolver(code,data);
					}
					else if(response.compare("GET RESULTS") == 0)
					{
						string message;
						response.clear();
						response = "OK! Please get ready to receive as BYTES Base64 Encoded!";
						sprintf(messageLength,"%8u",(int)response.size());
						message.append(messageLength);
						message.append(response);
						//printf("Sending: %s\n",message.c_str());
						//connection->sendLine(message);
						response.clear();
						//connection->recvLine(response);
						data = "RESULTS";
						code = 2;
						sendCommandToSolver(code,data);
					}
					else if(response.compare("STOP") == 0)
					{
						data = "STOP";
						code = 3;
						sendCommandToSolver(code,data);
					}
					else if(response.compare("BYE") == 0)
					{
						string message;
						response.clear();
						response = "Take care!";
						sprintf(messageLength,"%8u",(int)response.size());
						message.append(messageLength);
						message.append(response);
						connection->sendLine(message);
						connection->close();
						data = "STOP";
						code = 3;
						sendCommandToSolver(code,data);
						break;
					}
					else
					{
						string message;
						data += response ; data += " -- ";
						data += "Error! Command NOT recognized!";
						message.clear();sprintf(messageLength,"%8u",(int)data.size());
						message.append(messageLength);message.append(data);
						printf("Sending: %s\n",message.c_str());
						connection->sendLine(message);
						continue;
					}						
				}
				else
				{
					printf("Connection Lost!!\n");
					data = "STOP";
					code = 3;
					sendCommandToSolver(code,data);
					break;
				}
			}
		}
		catch(SocketException& )
		{
			printf("Connection Lost!!\n");
			string data = "STOP";
			int code = 3;
			sendCommandToSolver(code,data);
		}
		connection->close();
	}
}

string getRemoteCommand(int *code)
{
	string cmd("");
	if(pg.rank == 0)
	{
		#ifdef _WIN32
		WaitForSingleObject(remoteCommandsThreadMutexLock,INFINITE);
		#else
		pthread_mutex_lock(&remoteCommandsThreadMutexLock);
		#endif
		if(remoteCommands.size() > 0)
		{
			*code = remoteCodes[0];
			cmd.append(remoteCommands[0]);
			remoteCommands.erase(remoteCommands.begin());
			remoteCodes.erase(remoteCodes.begin());
		}
		#ifdef _WIN32
		ReleaseMutex(remoteCommandsThreadMutexLock);
		#else
		pthread_mutex_unlock(&remoteCommandsThreadMutexLock);
		#endif
	}
	return cmd;
}

   
void startRemoteCommandSocket()
{
	#ifdef _WIN32
	remoteCommandsThreadMutexLock = CreateMutex(NULL,FALSE,NULL);
	threadHd = CreateThread(0,0,(LPTHREAD_START_ROUTINE)listenToCommands,(LPVOID)NULL,0,&threadId);
	#else
	pthread_mutex_init(&remoteCommandsThreadMutexLock, NULL);
	pthread_create( &remoteCommandsThread, NULL, listenToCommands, (void*) NULL);
	#endif	    
	std::cout << "Remote Service Running....\r\nListening on Port: 30000\n";
}

void initMPI(int argc, char* argv[],int xNodes,int yNodes,bool remoteService = false)
{
	int myrank,size;
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	pg = getProcessGraph(myrank,xNodes,yNodes,size);
	
	if(myrank == 0 && remoteService)
	{
		startRemoteCommandSocket();
	}
}

void finMPI()
{
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
}

#endif
