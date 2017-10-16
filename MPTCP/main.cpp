#include <iostream>
#include <random>
#include <queue>
#include <vector>
#include <map>
#include <time.h>
#include <ctime>
#include <iomanip>
#include <map>
#include <random>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

class Packet;
int poissonFunction(int mean);
unsigned int numPackets=5;
int timeSimulation=1;
double timeSlot=0.0001;
long numberOfTimeSlot=(timeSimulation/timeSlot);
//unsigned short int maxNumberSubflow=5;
unsigned short int maxNumberNode=10;
unsigned short int maxNumberInterface;
unsigned short int maxNumberFlow=3;
unsigned long maxDataRate=1000000;
float probabilityInterface=0.5;
mt19937 initRandomEngine();
double generateUniform(mt19937 &engine);
class Node;
class Interface;
class NetworkQueue;
void assignInterfaceToNode(Interface &i,Node &n);
void sendPacket(Interface &interface,Node &node,long arrivalTimeSlot,vector<Node> vectorNodes);
class Flow;
void transportToSubflow(Node &node);
void transportToInterface(Node &node,vector<Interface> &interfaceVector);
void updateRoutingTable(Node &transmitter,Node &reciever,int cost,vector<Node> &vectorNodes);
void initilizeRoutingTable(Node &node,vector<Node> vectorNodes);
void showRoutingTable(vector<Node> vectorNodes);
unsigned int nodeNumber=1;
vector<vector<int>> routingTable;
void inistalRoutingTable( unsigned int nodeNumber);
void showMainRoutTbale();
void calculateBelmanFordAlgorithm(Node &node,vector<Node> &vectorNodes);
vector<vector<int>> distanceVectorAlgorithm(vector<Interface> interfaceVector,vector<Node> nodeVector);
vector<vector<int>> forwardingTable;
void printTwoDimentionalVector(vector<vector<int>> table);
class ApplicationLayer;
class TransportLayer;
class NetworkLayer;
class Subflow;
//void disableSubflow(Packet &packet,vector<Node> vectorNodes);
class Ack;


class Graph
{
        int V; // No. of vertices
        vector<int> *adj; // A dynamic array of adjacency lists
        void bridgeUtil(int v, bool visited[], int disc[], int low[],
                int parent[]);
    public:
        Graph(int V); // Constructor
        void addEdge(int v, int w); // function to add an edge to graph
        void bridge(); // prints all bridges
};

Graph::Graph(int V)
{
    this->V = V;
    adj = new vector<int> [V];
}

void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w);
    adj[w].push_back(v); // Note: the graph is undirected
}

void Graph::bridgeUtil(int u, bool visited[], int disc[], int low[],
        int parent[])
{
    // A static variable is used for simplicity, we can avoid use of static
    // variable by passing a pointer.
    static int time = 0;

    // Mark the current node as visited
    visited[u] = true;

    // Initialize discovery time and low value
    disc[u] = low[u] = ++time;

    // Go through all vertices aadjacent to this
    vector<int>::iterator i;
    for (i = adj[u].begin(); i != adj[u].end(); ++i)
    {
        int v = *i; // v is current adjacent of u

        // If v is not visited yet, then recur for it
        if (!visited[v])
        {
            parent[v] = u;
            bridgeUtil(v, visited, disc, low, parent);

            // Check if the subtree rooted with v has a connection to
            // one of the ancestors of u
            low[u] = min(low[u], low[v]);

            // If the lowest vertex reachable from subtree under v is
            // below u in DFS tree, then u-v is a bridge
            if (low[v] > disc[u])
                cout << u << " " << v << endl;
        }

        // Update low value of u for parent function calls.
        else if (v != parent[u])
            low[u] = min(low[u], disc[v]);
    }
}

// DFS based function to find all bridges. It uses recursive function bridgeUtil()
void Graph::bridge()
{
    // Mark all the vertices as not visited
    bool *visited = new bool[V];
    int *disc = new int[V];
    int *low = new int[V];
    int *parent = new int[V];

    // Initia4lize parent and visited arrays
    for (int i = 0; i < V; i++)
    {
        parent[i] = 0;
        visited[i] = false;
    }

    // Call the recursive helper function to find Bridges
    // in DFS tree rooted with vertex 'i'
    for (int i = 0; i < V; i++)
        if (visited[i] == false)
            bridgeUtil(i, visited, disc, low, parent);
}
class Packet{
public:

    unsigned long sequenceNumber;
    float percentSend;
    unsigned short int sourceAddress;
    unsigned short int destinationAddress;
    double arrivaltime;
    pair<unsigned int,unsigned int> packetFlow;
    unsigned short int subflowID;
    unsigned short int flowID;
    bool packetLoss=false;
    bool acknowledgment=false;
    Packet() {}

};
class ApplicationLayer{
public:

    queue<Packet> applicationQueue;
    void generatePacket(unsigned long numberPacket,Node &node,long timeSlotIndex);



};

class TransportLayer{
public:

    vector<Flow> flows;
    vector<AckFlow> ackFlows;
    unsigned short int numberFlows=0;
    void generateSubflow(Node &node);
    void appToTrnaspor(Node &node);
    void flowToSubflow(Node &node);
    void initSubflow(Flow &flow);
    void acknowledgement(Packet &packet);
    void generateAckFlow(Node &node,vector<Node> &vectorNodes);
    void handleAck();
    void sendAck(Node &node,Packet &packet);
};
class NetworkLayer{
public:

    vector<NetworkQueue> networkQueues;
    void transportToNetwork(Node &node);
    vector<NetworkQueue> recieveQueuesNetworkLayer;
    void generateRecievQueue(Interface &interface,Node &node);
    void detectPacketNetworkQueue(Packet &packet,Node &node);
};

class Flow{
public:
    queue<Packet> mainFlow;
    unsigned int sourecFlow;
    unsigned int destinationFlow;
    unsigned int numberSubflow;
    vector<Subflow> subflows;
    unsigned long flowSize;
    unsigned long lastPacket=0;
    bool isMultipath=true;
    unsigned short int id;
};
class AckFlow{
public:
    vector<Packet> mainFlow;
    unsigned int sourecFlow;
    unsigned int destinationFlow;
    unsigned int numberSubflow;
    vector<AckSubFlow> ackSubflow;
    unsigned long flowSize;
    unsigned long lastPacket=0;
    bool isMultipath=true;
    unsigned short int id;
};



class Node
{

    public:
    ///Vector of all Queues in the Network Layer of a Node
    vector<NetworkQueue> networkQueues;
    unsigned int id ;
    int positionX;
    int positionY;
    vector<Interface> connectedInterfaces;
    /// a Queue for all packets that their destination is this Node
    queue<Packet> reciveQueue;
    /// Flow that Creates in the Transport Layer
    Flow flow;
    ///the First item Node that Connected  second item interface
    vector<pair<int,int>> strightNeighbors;///First Node Second Interface
    vector<vector<int>> routingTable;
    //queue<Packet> forwardingQueue;
    ApplicationLayer appLayer;
    TransportLayer tcpLayer;
    NetworkLayer networkLayer;



};

class NetworkQueue{
public:
    int id;
    int interface;
    queue<Packet> q;
    Node node;
    pair<int,int> sourceDestQueue;
};

class Interface
{
public:

    pair<unsigned short int,unsigned short int> connectedNode;
    queue<Packet> interfaceQueue;
    int id;
    unsigned long dataRate;
    double numSendPacketTimeSlot;
    unsigned short int weight=1;
    NetworkQueue networkQueue;
};


class Subflow{
public:
    unsigned short int id;
    queue<Packet> q;
    Flow flow;
    unsigned int  sourceNode;
    unsigned int destinationNode;
    Interface interface;

};
class AckSubFlow{
public:
    unsigned short int id;
    vector<Packet> q;
    AckFlow ackFlow;
    unsigned int sourceNode;
    unsigned int destinationNode;
    Interface interface;
};

int randomIntNumber(unsigned int number){
    unsigned int randomNumber=rand()%number;
    return randomNumber;
}

void randomTopology(){
    nodeNumber=randomIntNumber(maxNumberNode);
    while(nodeNumber==0||nodeNumber==1){
        nodeNumber=randomIntNumber(maxNumberNode);
    }
    cout<<"Node number is : "<<nodeNumber<<"\n"<<"\n";
    vector<Node> vectorNodes;
    for(int i=0;i<nodeNumber;i++){
        Node node;
        node.id=i;
        node.positionX=randomIntNumber(100);
        node.positionY=randomIntNumber(100);
        cout<<"node "<<i<<" positionX: "<<node.positionX<<"\n";
        cout<<"node "<<i<<" positionY: "<<node.positionY<<"\n";
        vectorNodes.push_back(node);
        cout<<"\n";
    }
    ///Initilize the Flows
    for(int k=0;k<vectorNodes.size();k++){
        unsigned short temp=randomIntNumber(maxNumberFlow);
        if(temp==0){
            temp=1;
        }
        for(int i=0;i<temp;i++){
            Flow flow;
            flow.id=i;
            flow.sourecFlow=k;
            flow.destinationFlow=randomIntNumber(vectorNodes.size());
            if(k!=flow.destinationFlow){
                vectorNodes.at(k).tcpLayer.flows.push_back(flow);
                cout<<"Node "<<k<<"has a flow("<<k<<","<<flow.destinationFlow<<")"<<"\n";
            }
        }
        vectorNodes.at(k).tcpLayer.numberFlows=vectorNodes.at(k).tcpLayer.flows.size();
        cout<<"Node "<<k<<" has "<<vectorNodes.at(k).tcpLayer.numberFlows<<" Flows"<<"\n";

    }
    cout<<"Setup nodes is finished!"<<"\n";
    unsigned short int numInterface=0;
    Graph graph(nodeNumber);

    vector<Interface> interfaceVector;
    mt19937 re=initRandomEngine();
    for(int i=0;i<vectorNodes.size();i++){
        for(int j=i+1;j<vectorNodes.size();j++){

                if(generateUniform(re)>probabilityInterface){
                    Interface interface;
                    interface.id=numInterface;
                    numInterface++;
                    interface.connectedNode.first=i;
                    interface.connectedNode.second=j;
                    unsigned long tempDataRate=randomIntNumber(maxDataRate);
                    interface.dataRate=tempDataRate;
                    interface.numSendPacketTimeSlot=interface.dataRate*timeSlot;
                    cout<<"Interface "<<interface.id<<" connected to node: "<<vectorNodes.at(i).id;
                    cout<<" and node: "<<vectorNodes.at(j).id;
                    cout<<" DataRate : "<<interface.dataRate<<"\n";
                    interfaceVector.push_back(interface);
                    vectorNodes.at(i).connectedInterfaces.push_back(interface);
                    pair<int,int> temp;
                    temp.first=j;
                    temp.second=interface.id;
                    vectorNodes.at(i).strightNeighbors.push_back(temp);
                    NetworkQueue q;
                    q.interface=interface.id;
                    q.node=vectorNodes.at(i);
                    q.sourceDestQueue.first=i;
                    q.sourceDestQueue.second=j;
                    vectorNodes.at(i).networkQueues.push_back(q);
                    interface.networkQueue=q;
                    graph.addEdge(i,j);
                    vectorNodes.at(j).networkLayer.generateRecievQueue(interface,vectorNodes.at(j));




                    ///Full Duplexing
                    Interface interface1;
                    interface1.id=numInterface;
                    numInterface++;
                    interface1.connectedNode.first=j;
                    interface1.connectedNode.second=i;
                    interface1.dataRate=tempDataRate;
                    interface1.numSendPacketTimeSlot=interface1.dataRate*timeSlot;
                    cout<<"Interface "<<interface1.id<<" connected to node: "<<vectorNodes.at(j).id;
                    cout<<" and node: "<<vectorNodes.at(i).id;
                    cout<<" DataRate : "<<interface1.dataRate<<"\n";
                    interfaceVector.push_back(interface1);
                    vectorNodes.at(j).connectedInterfaces.push_back(interface1);
                    pair<int,int> temp1;
                    temp1.first=i;
                    temp1.second=interface1.id;
                    vectorNodes.at(j).strightNeighbors.push_back(temp1);
                    NetworkQueue q1;
                    q1.interface=interface1.id;
                    q1.node=vectorNodes.at(j);
                    q1.sourceDestQueue.first=j;
                    q1.sourceDestQueue.second=i;
                    vectorNodes.at(j).networkQueues.push_back(q1);
                    interface1.networkQueue=q1;
                    graph.addEdge(j,i);
                    vectorNodes.at(i).networkLayer.generateRecievQueue(interface1,vectorNodes.at(i));


                }

        }

    }


    cout<<"-------------------------------------"<<"\n";
    if(numInterface==0){
        cout<<"We have no interface !"<<"\n";
    }
    ///init the subflows per Node
    for(int i=0;i<vectorNodes.size();i++){
        if(!vectorNodes.at(i).connectedInterfaces.empty()){
            vectorNodes.at(i).tcpLayer.generateSubflow(vectorNodes.at(i));
        }

    }




    forwardingTable=distanceVectorAlgorithm(interfaceVector,vectorNodes);
    cout<<"Printing Forwarding Table"<<"\n";
    printTwoDimentionalVector(forwardingTable);
    cout<<"Graph Bridge is :"<<"\n";
    graph.bridge();


   ///Transferring Packets in each Time Slot
    for(long r=0;r<1;r++){
        ///Application Layer Generate the Packets
        for(int t=0;t<vectorNodes.size();t++){
            vectorNodes.at(t).appLayer.generatePacket(numPackets,vectorNodes.at(t),r);
        }
        ///transferring packets from Application Layer to the Transport Layer flows
        for(int u=0;u<vectorNodes.size();u++){
            if(!vectorNodes.at(u).appLayer.applicationQueue.empty()){
                vectorNodes.at(u).tcpLayer.appToTrnaspor(vectorNodes.at(u));
            }


        }

        ///Transferring Packets from flows to the Subflows;
        for(int n=0;n<vectorNodes.size();n++){
            if(vectorNodes.at(n).tcpLayer.numberFlows!=0){
                vectorNodes.at(n).tcpLayer.flowToSubflow(vectorNodes.at(n));
            }

        }

        ///Tansport the Packets to the Network Layer
        for(int p=0;p<vectorNodes.size();p++){
            vectorNodes.at(p).networkLayer.transportToNetwork(vectorNodes.at(p));
        }



        vector<double> intfcTime;
        for(int g=0;g<interfaceVector.size();g++){
            double f=interfaceVector.at(g).numSendPacketTimeSlot;
            intfcTime.push_back(f);
        }
        ///Transmit packets to interface Queue
        for(int l=0;l<vectorNodes.size();l++){
            transportToInterface(vectorNodes.at(l),interfaceVector);
        }

        for(int w=0;w<interfaceVector.size();w++){
        cout<<"Interface "<<w<<"has a queue size of : "<<interfaceVector.at(w).interfaceQueue.size()<<"\n";
        }

        for(int h=0;h<intfcTime.size();h++){
            unsigned short int recieverIndex=interfaceVector.at(h).connectedNode.second;
            while(intfcTime.at(h)>0){
                if(intfcTime.at(h)>=1){
                    if(!interfaceVector.at(h).interfaceQueue.empty()){
                        if(interfaceVector.at(h).interfaceQueue.front().percentSend==0){
                            sendPacket(interfaceVector.at(h),vectorNodes.at(recieverIndex),r,vectorNodes);
                            intfcTime.at(h)=intfcTime.at(h)-1;
                            }
                        else{
                            int temp=interfaceVector.at(h).interfaceQueue.front().percentSend;
                            sendPacket(interfaceVector.at(h),vectorNodes.at(recieverIndex),r,vectorNodes);
                            interfaceVector.at(h).interfaceQueue.front().percentSend=temp;
                            intfcTime.at(h)=intfcTime.at(h)-1;
                        }


                    }
                    else{
                        intfcTime.at(h)=intfcTime.at(h)-1;
                    }

                }
                else{
                    if(!interfaceVector.at(h).interfaceQueue.empty()){
                        interfaceVector.at(h).interfaceQueue.front().percentSend+=intfcTime.at(h)*100;
                        if(interfaceVector.at(h).interfaceQueue.front().percentSend>=100){
                            int temp=interfaceVector.at(h).interfaceQueue.front().percentSend-100;
                            sendPacket(interfaceVector.at(h),vectorNodes.at(recieverIndex),r,vectorNodes);
                            interfaceVector.at(h).interfaceQueue.front().percentSend+=temp;
                            intfcTime.at(h)=intfcTime.at(h)-1;
                        }
                        else{
                            intfcTime.at(h)=intfcTime.at(h)-1;
                        }

                    }
                    else{
                        intfcTime.at(h)=intfcTime.at(h)-1;
                    }
                }
            }

        }
        for(int b=0;b<vectorNodes.size();b++){
        cout<<"Recive Queue in Node: "<<b<<" is: "<<vectorNodes.at(b).reciveQueue.size()<<"\n";
    }

    }

    cout << "Following are connected components \n";
//    graph.connectedComponents();
}









mt19937 initRandomEngine(){
    std::random_device rd;
    std::mt19937 e2(rd());
    e2.seed(time(NULL));
    return e2;
}


double generateUniform(mt19937 &engine){

    uniform_real_distribution<> dist(0, 1);
    return dist(engine);
}·
void sendPacket(Interface &interface,Node &node,long arriveTimeSlot,vector<Node> vectorNodes){
    Packet packet=interface.interfaceQueue.front();
    interface.interfaceQueue.pop();
    if(packet.destinationAddress==node.id){
        packet.percentSend=0;
        packet.arrivaltime=timeSlot*(arriveTimeSlot-packet.arrivaltime);
        for(int i=0;i<node.networkLayer.recieveQueuesNetworkLayer.size();i++){
            if(node.networkLayer.recieveQueuesNetworkLayer.at(i).interface==interface.id){
                node.networkLayer.recieveQueuesNetworkLayer.at(i).q.push(packet);
                node.tcpLayer.handleAck();
            }
        }
    }
    else{
        Packet packet;
        packet=interface.interfaceQueue.front();
        interface.interfaceQueue.pop();
        int source=node.id;
        int destination=packet.destinationAddress;
        int forwarded=forwardingTable.at(source).at(destination);
        if(forwarded==-1){
            ///Packet Drop and tell the Sender
        }
        else{
            for(int i=0;i<node.networkQueues.size();i++){
                if(node.networkQueues.at(i).sourceDestQueue.first==node.id&&
                    node.networkQueues.at(i).sourceDestQueue.second==destination){
                    node.networkQueues.at(i).q.push(packet);
                }
            }
        }
    }

}

void transportToInterface(Node &node,vector<Interface> &interfaceVector){
    if(!node.networkQueues.empty()){
        for(int i=0;i<node.networkQueues.size();i++){
            while(!node.networkQueues.at(i).q.empty()){
                interfaceVector.at(node.networkQueues.at(i).interface).interfaceQueue.push(node.networkQueues.at(i).q.front());
                node.networkQueues.at(i).q.pop();
            }

        }
    }

}
vector<vector<int>> distanceVectorAlgorithm(vector<Interface> interfaceVector,vector<Node> nodeVector){
    int graph[nodeNumber][nodeNumber];
    int i,j,k,t=1;
    int nn=nodeNumber;

    /* Initialize graph*/
    for (i=0;i<nn;i++)
    {
        for(j=0;j<nn;j++)
        {
            graph[i][j]=-1;
        }
    }

    /* Vertex names */

    char ch1[26]={'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W',
    'X','Y','Z'};
    char ch[nodeNumber];
    for(int i=0;i<nodeNumber;i++){
        ch[i]=ch1[i];
    }
    /*char ch[26]={'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W',
    'X','Y','Z'};*/


    /* Get input */
    for (i=0;i<nn;i++)
    {
        for(j=0;j<nn;j++)
        {
            if(i==j)
            {
                graph[i][j]=0;
            }
            if(graph[i][j]==-1)
            {
                if(!nodeVector.at(i).strightNeighbors.empty()){
                    for(int z=0;z<nodeVector.at(i).strightNeighbors.size();z++){
                        if(nodeVector.at(i).strightNeighbors.at(z).first==j){
                            graph[i][j]=graph[j][i]=1;

                        }



                    }

                }

            }
        }
    }


    /* Initializing via */
    int via[nodeNumber][nodeNumber];
    for (i=0;i<nn;i++)
    {
        for(j=0;j<nn;j++)
        {
            via[i][j]=-1;
        }
    }

    cout<<"\n After Initialization";
    /* Display table initialization */
    for (i=0;i<nn;i++)
    {
        cout<<"\n"<<ch[i]<<" Table";
        cout<<"\nNode\tDist\tVia";
        for(j=0;j<nn;j++)
        {
            cout<<"\n"<<ch[j]<<"\t"<<graph[i][j]<<"\t"<<via[i][j];
        }
    }

    //sharing table
    int sh[50][50][50];
    for(i=0;i<nn;i++)
    {
        for(j=0;j<nn;j++)
        {
            for (k=0;k<nn;k++)
            {
                /* Check if edge is present or not*/
                if((graph[i][j]>-1)&&(graph[j][k]>-1))
                {
                    sh[i][j][k]=graph[j][k]+graph[i][j];
                }
                else
                {
                    sh[i][j][k]=-1;
                }
            }
        }
    }

    /*displaying shared table */
    for(i=0;i<nn;i++)
    {
        cout<<"\n\n For "<<ch[i];
        for (j=0;j<nn;j++)
        {
            cout<<"\n From "<<ch[j];
            for(k=0;k<nn;k++)
            {
                cout<<"\n "<<ch[k]<<" "<<sh[i][j][k];
            }
        }
    }

    /* Updating */
    int final[50][50];
    for(i=0;i<nn;i++)
    {
        for(j=0;j<nn;j++)
        {
            /* Copy initial value from input graph*/
            final[i][j]=graph[i][j];
            via[i][j]=j;


            /*Update them*/
            /* Check condition a - b - c */
            for(k=0;k<nn;k++)
            {

                if((final[i][j]>sh[i][k][j]) || (final[i][j] == -1))
                {
                    if(sh[i][k][j]>-1)
                    {
                        final[i][j]=sh[i][k][j];
                        via[i][j]=k;
                    }
                }
            }
            /* After considering three vertex if final not found
                consider 4th
                a- b- c- d
            */

            if(final[i][j]==-1)
            {
                for(k=0;k<nn;k++)
                {

                    if((final[i][k]!=-1)&&(final[k][j]!=-1))
                    {
                        if((final[i][j]==-1) || ((final[i][j]!=-1) &&(final[i][j]>final[i][k]+final[k][j])))
                        {
                            if(final[i][k]+final[k][j]>-1)
                            {
                                final[i][j]=final[i][k]+final[k][j];
                                via[i][j]=-1;
                            }
                        }
                    }

                }
            }
        }
    }

    cout<<"\n After Update :";
    /* Display table Updation */
    for (i=0;i<nn;i++)
    {
        cout<<"\n"<<ch[i]<<" Table";
        cout<<"\nNode\tDist\tVia";
        for(j=0;j<nn;j++)
        {
            cout<<"\n"<<ch[j]<<"\t"<<final[i][j]<<"\t";
            if(i==via[i][j])
                cout<<"-";
            else
                cout<<ch[via[i][j]];
        }
    }
    cout<<"\n";
    cout<<"-----------------------------------"<<"\n";
    vector<vector<int>> temp;
    for(int e=0;e<nodeNumber;e++){
        vector<int> rowTemp;
        for(int t=0;t<nodeNumber;t++){
            rowTemp.push_back(via[e][t]);
            if(e==t){
                rowTemp.at(t)=e;
            }
        }
        temp.push_back(rowTemp);
    }
    return temp;
}
int poissonFunction(int mean){
    default_random_engine generator;
    poisson_distribution<int> distribution(mean);
    return distribution(generator);
}
void printTwoDimentionalVector(vector<vector<int>> table){
    for(int i=0;i<forwardingTable.size();i++){
        for(int j=0;j<forwardingTable.at(i).size();j++){
            cout<<forwardingTable.at(i).at(j)<<" ";
        }
        cout<<"\n";
    }
    cout<<"\n";

}

void ApplicationLayer::generatePacket(unsigned long numberPacket,Node &node,long timeSlotIndex){
    if(!node.tcpLayer.flows.empty()){
        for(int i=0; i<numberPacket; i++){
                Packet packet;
                packet.sourceAddress=node.id;
                packet.percentSend=0;
                packet.arrivaltime=timeSlotIndex;
                applicationQueue.push(packet);
        }
        cout<<"Node "<<node.id<<" in the Application Layer has : "<<applicationQueue.size()<<" Packets"<<"\n";
    }
}
void TransportLayer::generateSubflow(Node &node,vector<Node> vectorNodes){
    if(node.tcpLayer.numberFlows!=0&&node.connectedInterfaces.size()!=0){
        for(int i=0;i<node.tcpLayer.numberFlows;i++){
            if(flows.at(i).isMultipath){
                for(int j=0;j<node.connectedInterfaces.size();j++){
                    Subflow subflow;
                    subflow.destinationNode=flows.at(i).destinationFlow;
                    subflow.flow=flows.at(i);
                    subflow.id=j;
                    subflow.interface=node.connectedInterfaces.at(j);
                    subflow.sourceNode=node.id;
                    flows.at(i).subflows.push_back(subflow);
                    flows.at(i).numberSubflow++;


                }
            }
        }
    }

}
void TransportLayer::appToTrnaspor(Node &node){
    if(!node.tcpLayer.flows.empty()){
        int appSize=node.appLayer.applicationQueue.size();
        int temp=0;
        if(appSize!=0){
            while(appSize>0){
                for(int i=0;i<node.tcpLayer.numberFlows;i++){
                    if(appSize>0){
                            Packet packet=node.appLayer.applicationQueue.front();
                            node.appLayer.applicationQueue.pop();
                            packet.destinationAddress=node.tcpLayer.flows.at(i).destinationFlow;
                            packet.packetFlow.first=node.id;
                            packet.packetFlow.second=node.tcpLayer.flows.at(i).destinationFlow;
                            packet.percentSend=0;
                            packet.sequenceNumber=node.tcpLayer.flows.at(i).lastPacket;
                            node.tcpLayer.flows.at(i).lastPacket++;
                            packet.sourceAddress=node.id;
                            node.tcpLayer.flows.at(i).mainFlow.push(packet);
                            appSize--;
                    }

                }
            }

        }
    }





    if(node.tcpLayer.numberFlows!=0){
        for(int z=0;z<node.tcpLayer.numberFlows;z++){
            node.tcpLayer.flows.at(z).flowSize=node.tcpLayer.flows.at(z).mainFlow.size();
        }
        for(int i=0;i<node.tcpLayer.flows.size();i++){
          cout<<"Node "<<node.id<<" in the flow "<<i<<" has "<<flows.at(i).mainFlow.size()<<" Packets"<<"\n";
        }
    }




}
void TransportLayer::flowToSubflow(Node &node){
    if(!node.tcpLayer.flows.empty()){
        for(int j=0;j<node.tcpLayer.flows.size();j++){
            unsigned long sizeFlow=node.tcpLayer.flows.at(j).mainFlow.size();
            if(!node.tcpLayer.flows.at(j).subflows.empty()){
                while(sizeFlow>0){
                    for(int i=0;i<node.tcpLayer.flows.at(j).subflows.size();i++){
                        if(sizeFlow>0){
                            Packet packet=node.tcpLayer.flows.at(j).mainFlow.front();
                            node.tcpLayer.flows.at(j).mainFlow.pop();
                            packet.subflowID=i;
                            packet.flowID=node.tcpLayer.flows.at(j).id;
                            node.tcpLayer.flows.at(j).subflows.at(i).q.push(node.tcpLayer.flows.at(j).mainFlow.front());
                            sizeFlow--;
                        }
                        cout<<"Node "<<node.id<<" in flow "<<j<<" in Subflow "<<i<<"has "<<node.tcpLayer.flows.at(j).subflows.at(i).q.size()<<" Packets"<<"\n";
                    }


                }
            }
        }
    }


}
void NetworkLayer::transportToNetwork(Node &node){
    for(int i=0;i<node.tcpLayer.numberFlows;i++){
        unsigned long temp=0;
        while(temp<node.tcpLayer.flows.at(i).flowSize){
            for(int j=0;j<node.tcpLayer.flows.at(i).subflows.size();j++){
                if(!node.tcpLayer.flows.at(i).subflows.at(j).q.empty()){
                    Packet packet=node.tcpLayer.flows.at(i).subflows.at(j).q.front();
                    node.tcpLayer.flows.at(i).subflows.at(j).q.pop();
                    //unsigned short int destination=packet.destinationAddress;
                    //unsigned short int source=packet.sourceAddress;
                    // int forwarded=forwardingTable.at(source).at(destination);
                    for(int w=0;w<node.connectedInterfaces.size();w++){
                        if(node.tcpLayer.flows.at(i).subflows.at(j).interface.id==node.connectedInterfaces.at(w).id){
                            for(int z=0;z<node.networkQueues.size();z++){
                                if(node.networkQueues.at(z).interface==node.connectedInterfaces.at(w).id){
                                    node.networkQueues.at(z).q.push(packet);
                                }
                            }
                        }

                      /*if(forwarded==node.connectedInterfaces.at(w).connectedNode.second){
                        for(int z=0;z<node.networkQueues.size();z++){
                            if(node.networkQueues.at(z).interface==node.connectedInterfaces.at(w).id){
                                node.networkQueues.at(z).q.push(packet);
                            }
                        }
                      }*/
                    }
                }
            }
            temp++;
        }

   }
   for(int i=0;i<node.networkQueues.size();i++){
    cout<<"Node "<<node.id<<" in network queue "<<i<<" has : "<<node.networkQueues.at(i).q.size()<<" Packets"<<"\n";
   }

}

void TransportLayer::generateAckFlow(Node &node,vector<Node> vectorNodes){
    for(int i=0;i<node.tcpLayer.flows.size();i++){
        AckFlow ackflow;
        ackflow.id=node.tcpLayer.flows.at(i).id;
        ackflow.isMultipath=node.tcpLayer.flows.at(i).isMultipath;
        ackflow.numberSubflow=node.tcpLayer.flows.at(i).numberSubflow;
        ackflow.sourecFlow=node.tcpLayer.flows.at(i).sourecFlow;
        ackflow.destinationFlow=node.tcpLayer.flows.at(i).destinationFlow;
        for(int j=0;j<ackflow.numberSubflow;j++){
            AckSubFlow ackSubflow;
            ackSubflow.id=j;
            ackSubflow.destinationNode=ackflow.destinationFlow;
            ackSubflow.ackFlow=ackflow;
            ackSubflow.sourceNode=node.id;
            ackflow.ackSubflow.push_back(ackSubflow);
        }
        unsigned short int reciever=ackflow.destinationFlow;
        vectorNodes.at(reciever).tcpLayer.ackFlows.push_back(ackflow);
        cout<<"The Ack Flow in Node "<<reciever<<" with "<<ackflow.numberSubflow<<" Ack Subflow initial"<<"\n";
    }
}
void NetworkLayer::generateRecievQueue(Interface &interface,Node &node){
    NetworkQueue recieveQueue;
    recieveQueue.id=interface.id;
    recieveQueue.interface=interface.id;
    recieveQueue.node=node;
    node.networkLayer.recieveQueuesNetworkLayer.push_back(recieveQueue);
}
void TransportLayer::handleAck(Node &node){
    for(int i=0;i<node.networkLayer.recieveQueuesNetworkLayer.size();i++){
        while(!node.networkLayer.recieveQueuesNetworkLayer.at(i).q.empty()){
            Packet packet=node.networkLayer.recieveQueuesNetworkLayer.at(i).q.front();
            node.networkLayer.recieveQueuesNetworkLayer.at(i).q.pop();
            for(int j=0;j<node.tcpLayer.ackFlows.size();j++){
                if(packet.flowID==node.tcpLayer.ackFlows.at(j).id){
                    for(int z=0;z<node.tcpLayer.ackFlows.at(j).ackSubflow.size();z++){
                        if(packet.subflowID==node.tcpLayer.ackFlows.at(j).ackSubflow.at(z).id){
                            int lastElement=node.tcpLayer.ackFlows.at(j).ackSubflow.at(z).q.size();
                            if(packet.sequenceNumber!=lastElement){
                                ///Packet Loss
                                cout<<"PacketLoss Detected on Packet With Sequence Number: ";
                                cout<<packet.sequenceNumber<<"\n";
                                node.tcpLayer.packetLoss(node,packet);
                            }
                            else{
                                node.tcpLayer.ackFlows.at(j).ackSubflow.at(z).q.push_back(packet);
                                ///Send The Acknowledgment
                                node.tcpLayer.sendAck(node,packet);
                            }
                        }
                    }
                }
            }
        }
    }
}
void TransportLayer::packetLoss(Node &node,Packet &packet){
    unsigned short int sourceAddress=packet.sourceAddress;
    unsigned short int destinationAddress=packet.destinationAddress;
    packet.destinationAddress=sourceAddress;
    packet.sourceAddress=destinationAddress;
    packet.packetLoss=true;
    packet.percentSend=0;
    ///Send the Packet to the Network Queue
    node.networkLayer.detectPacketNetworkQueue(node,packet);
}
void NetworkLayer::detectPacketNetworkQueue(Packet &packet,Node &node){
    unsigned int forwardedTo=forwardingTable.at(packet.sourceAddress).at(packet.destinationAddress);
    for(int i=0;i<node.connectedInterfaces.size();i++){
        if(node.connectedInterfaces.at(i).connectedNode.first==node.id
           &&node.connectedInterfaces.at(i).connectedNode.second==forwardedTo){
            for(int j=0;node.networkLayer.networkQueues.size();j++){
                if(node.connectedInterfaces.at(i).id==node.networkLayer.networkQueues.at(j).interface){
                    node.networkLayer.networkQueues.at(j).q.push(packet);
                }
            }
        }
    }
}
void TransportLayer::sendAck(Node &node,Packet &packet){
    packet.acknowledgment=true;
    unsigned short int source=packet.destinationAddress;
    unsigned short int destination=packet.source;
    packet.destinationAddress=destination;
    packet.sourceAddress=source;
    packet.percentSend=0;
    ///Send To Network Layer
    node.networkLayer.detectPacketNetworkQueue(packet,node);
    cout<<"Ack is Sended "<<"\n";
}

int main()
{

    srand(time(NULL));
    randomTopology();
    return 0;
}

