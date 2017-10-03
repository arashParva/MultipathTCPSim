#include <iostream>
#include <random>
#include <queue>
#include <vector>
#include <map>
#include <time.h>
#include <ctime>
using namespace std;
class Interface;
class Packet;
long numPackets=10000;
int timeSimulation=1;
double timeSlot=0.0001;
long numberOfTimeSlot=(timeSimulation/timeSlot);
unsigned short int maxNumberSubflow=10;
unsigned short int maxNumberNode=10;
unsigned short int maxNumberInterface;

class Node
{

    public:
    queue<Packet> networkQueue;
    int num_queues;
    std::vector<queue<Packet>>  nodeQueues;
    long unsigned int id ;
    int positionX;
    int positionY;

};
class Packet
{
public:

    unsigned int sequenceNumber;
    unsigned short int sendPercent;
    float percentSend;
    Packet() {}

};

class Interface
{
public:

    queue<Packet> interfaceQueue;
    long unsigned int id;
    double delay;
    unsigned long dataRate;
    double numSendPacketTimeSlot;

};


void generatePacket(int numberPacket,queue<Packet> &nodeQueue)
{
    for(int i=0; i<numberPacket; i++)
    {
        Packet packet;
        packet.sendPercent=0;
        nodeQueue.push(packet);
    }

}
void simpleTopology()
{
    ///define Nodes;
    Node node1;
    Node node2;
    Node node3;
    Node node4;


    ///initialize number of subflows;
    node1.num_queues=2;
    node2.num_queues=1;
    node3.num_queues=1;
    node4.num_queues=1;

    ///Initialize the Queue nodes
    for(int i=0; i<node1.num_queues; i++)
    {
        queue<Packet> nodeQueue;
        node1.nodeQueues.push_back(nodeQueue);
    }
    for(int j=0; j<node2.num_queues; j++)
    {
        queue<Packet> nodeQueue;
        node2.nodeQueues.push_back(nodeQueue);
    }
    for(int z=0; z<node3.num_queues; z++)
    {
        queue<Packet> nodeQueue;
        node3.nodeQueues.push_back(nodeQueue);
    }
    for(int f=0;f<node4.num_queues;f++){
        queue<Packet> nodeQueue;
        node4.nodeQueues.push_back(nodeQueue);
    }

    ///create interfaces
    Interface interface1;
    Interface interface2;
    Interface interface3;

    ///initial DataRate Interface
    interface1.dataRate=4096;
    interface2.dataRate=2048;
    interface3.dataRate=1024;

    interface1.numSendPacketTimeSlot=interface1.dataRate*timeSlot;
    interface2.numSendPacketTimeSlot=interface2.dataRate*timeSlot;
    interface3.numSendPacketTimeSlot=interface3.dataRate*timeSlot;

    ///generate packets
    generatePacket(numPackets,node1.nodeQueues.at(0));
    generatePacket(numPackets,node1.nodeQueues.at(1));
    generatePacket(numPackets,node2.nodeQueues.at(0));

    ///deliver packets from transport queue to network queue
    while(!node1.nodeQueues.at(0).empty()&& !node1.nodeQueues.at(1).empty())
    {
        node1.networkQueue.push(node1.nodeQueues.at(0).front());
        node1.nodeQueues.at(0).pop();
        node1.networkQueue.push(node1.nodeQueues.at(1).front());
        node1.nodeQueues.at(1).pop();
    }

    while(!node2.nodeQueues.at(0).empty()){

        node2.networkQueue.push(node2.nodeQueues.at(0).front());
        node2.nodeQueues.at(0).pop();
    }
    float i1;
    float i2;
    float i3;
    ///transferring packets
    for(long z=0;z<numberOfTimeSlot;z++){
        i1=interface1.numSendPacketTimeSlot;
        while(i1>0){
            if(i1>=1){
                if(!node1.networkQueue.empty()){
                            interface1.interfaceQueue.push(node1.networkQueue.front());
                            node1.networkQueue.pop();
                            node3.networkQueue.push(interface1.interfaceQueue.front());
                            cout<<"SendPacket"<<"\n";
                            interface1.interfaceQueue.pop();
                            i1=i1-1;
                }

            }
            else{
                if(!node1.networkQueue.empty()){
                            float temp=i1*100;
                            node1.networkQueue.front().percentSend+=temp;
                            i1=i1-1;
                            if(node1.networkQueue.front().percentSend>=100){
                                interface1.interfaceQueue.push(node1.networkQueue.front());
                                node1.networkQueue.pop();
                                node3.networkQueue.push(interface1.interfaceQueue.front());
                                cout<<"SendPacket"<<"\n";
                                interface1.interfaceQueue.pop();
                                node3.networkQueue.back().percentSend=0;
                            }
                }
            }
        }
        i2=interface2.numSendPacketTimeSlot;
        while(i2>0){
            if(i2>=1){
                if(!node2.networkQueue.empty()){

                            interface2.interfaceQueue.push(node2.networkQueue.front());
                            node2.networkQueue.pop();
                            node3.networkQueue.push(interface2.interfaceQueue.front());
                            cout<<"SendPacket"<<"\n";
                            interface2.interfaceQueue.pop();
                            i2=i2-1;
                }

            }
            else{
                if(!node2.networkQueue.empty()){
                            node2.networkQueue.front().percentSend+=i2*100;
                            i2=i2-1;
                            if(node2.networkQueue.front().percentSend>=100){
                                node2.networkQueue.front().percentSend=0;
                                interface2.interfaceQueue.push(node2.networkQueue.front());
                                node2.networkQueue.pop();
                                node3.networkQueue.push(interface2.interfaceQueue.front());
                                cout<<"SendPacket"<<"\n";
                                interface2.interfaceQueue.pop();

                            }
                }

           }
    }




        i3=interface3.numSendPacketTimeSlot;
        while(i3>0){
            if(i3>=1){
                if(!node3.networkQueue.empty()){
                            interface3.interfaceQueue.push(node3.networkQueue.front());
                            node3.networkQueue.pop();
                            node4.networkQueue.push(interface3.interfaceQueue.front());
                            cout<<"SendPacket"<<"\n";
                            interface3.interfaceQueue.pop();
                            i3=i3-1;
                }
                i3=i3-1;
            }
            else{
                if(!node3.networkQueue.empty()){
                        if(node3.networkQueue.front().percentSend!=100){
                            node3.networkQueue.front().percentSend+=i3*100;
                            if(node3.networkQueue.front().percentSend>=100){
                                interface3.interfaceQueue.push(node3.networkQueue.front());
                                node3.networkQueue.pop();
                                node4.networkQueue.push(interface3.interfaceQueue.front());
                                cout<<"SendPacket"<<"\n";
                                node4.networkQueue.back().percentSend=0;
                                interface3.interfaceQueue.pop();
                            }
                            i3=i3-1;
                        }
                        else{
                            interface3.interfaceQueue.push(node3.networkQueue.front());
                            node3.networkQueue.pop();
                            node4.networkQueue.push(interface3.interfaceQueue.front());
                            cout<<"SendPacket"<<"\n";
                            interface3.interfaceQueue.pop();
                            node4.networkQueue.back().percentSend=0;
                            i3=i3-1;
                        }
                }

            i3=i3-1;
        }
    }
}

    cout<<"Number of packets in Node1 Queue: "<<node1.networkQueue.size()<<"\n";
    cout<<"Number of packets in Node2 Queue: "<<node2.networkQueue.size()<<"\n";
    cout<<"Number of packets in Node3 Queue: "<<node3.networkQueue.size()<<"\n";
    cout<<"Number of packets in Node4 Queue: "<<node4.networkQueue.size()<<"\n";
}
int randomIntNumber(unsigned int number){
    unsigned int randomNumber=rand()%number;
    return randomNumber;
}
void randomTopology(){
    unsigned int nodeNumber=randomIntNumber(maxNumberNode);
    while(nodeNumber==0){
        nodeNumber=randomIntNumber(maxNumberNode);
    }
    cout<<"Node number is : "<<nodeNumber<<"\n";
    vector<Node> vectorNodes;
    for(int i=0;i<nodeNumber;i++){
        Node node;
        node.id=i;
        node.positionX=randomIntNumber(100);
        node.positionY=randomIntNumber(100);
        cout<<"node "<<i<<" positionX: "<<node.positionX<<"\n";
        cout<<"node "<<i<<" positionY: "<<node.positionY<<"\n";
        unsigned short int numSubFlow=randomIntNumber(maxNumberSubflow);
        while(numSubFlow==0){
            numSubFlow=randomIntNumber(maxNumberSubflow);
        }
        cout<<"numSubflow node "<<i<<": is : "<<numSubFlow<<"\n";
        for(int j=0;j<numSubFlow;j++){
            queue<Packet> subFlowQueue;
            node.nodeQueues.push_back(subFlowQueue);
        }
        vectorNodes.push_back(node);

    }
    cout<<"Setup nodes is finished!"<<"\n";
    maxNumberInterface=nodeNumber*(nodeNumber-1);
    unsigned short int numInterfaces=randomIntNumber(maxNumberInterface);
    while(numInterfaces==0){
        numInterfaces=randomIntNumber(maxNumberInterface);
    }
    cout<<"Number of Interfaces is : "<<numInterfaces<<"\n";
    for(int m=0;m<numInterfaces;m++){
            Interface interface;
            interface.dataRate=1024;
            interface.id=i;
            interface.

    }




}

int main()
{
    srand (time(NULL));
    randomTopology();
    //simpleTopology();
    return 0;
}
