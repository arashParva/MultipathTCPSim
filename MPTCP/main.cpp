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

using namespace std;

class Packet;
int poissonFunction(int mean);
unsigned int numPackets=poissonFunction(5);
int timeSimulation=1;
double timeSlot=0.0001;
long numberOfTimeSlot=(timeSimulation/timeSlot);
unsigned short int maxNumberSubflow=5;
unsigned short int maxNumberNode=10;
unsigned short int maxNumberInterface;

unsigned long maxDataRate=1000000;
float probabilityInterface=0.5;
mt19937 initRandomEngine();
double generateUniform(mt19937 &engine);
class Node;
class Interface;
class NetworkQueue;
void assignInterfaceToNode(Interface &i,Node &n);
void sendPacket(Interface &interface,Node &node,long arrivalTimeSlot);
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

    // Initialize parent and visited arrays
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

class Flow{
public:
    queue<Packet> mainFlow;
    pair<unsigned int,unsigned int> transDestFlow;
    unsigned int numberSubflow;
    vector<queue<Packet>> subflows;
    unsigned long flowSize;
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

};
class Packet
{
public:

    unsigned int sequenceNumber;
    float percentSend;
    unsigned short int sourceAddress;
    unsigned short int destinationAddress;
    double arrivaltime;
    Packet() {}

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



void generatePacket(unsigned int numberPacket,Node &node,long timeSlotIndex)
{
    for(int i=0; i<numberPacket; i++)
    {
        Packet packet;
        packet.sequenceNumber=i;
        packet.sourceAddress=node.flow.transDestFlow.first;
        packet.destinationAddress=node.flow.transDestFlow.second;
        packet.percentSend=0;
        packet.arrivaltime=timeSlotIndex;
        node.flow.mainFlow.push(packet);
    }
    cout<<"Node "<<node.id<<" in main flow has "<<node.flow.mainFlow.size()<<" Packets"<<"\n";

}
int randomIntNumber(unsigned int number){
    unsigned int randomNumber=rand()%number;
    return randomNumber;
}
void transportToNetwork(Node &node){
    unsigned long temp=0;
    while(temp<node.flow.flowSize){
        for(int j=0;j<node.flow.numberSubflow;j++){
            if(!node.flow.subflows.at(j).empty()){
                Packet packet;
                packet=node.flow.subflows.at(j).front();
                node.flow.subflows.at(j).pop();
                unsigned short int destination=packet.destinationAddress;
                unsigned short int source=packet.sourceAddress;
                int forwarded=forwardingTable.at(source).at(destination);
                for(int i=0;i<node.connectedInterfaces.size();i++){
                    if(forwarded==node.connectedInterfaces.at(i).connectedNode.second){
                        for(int z=0;z<node.networkQueues.size();z++){
                            if(node.networkQueues.at(z).interface==node.connectedInterfaces.at(i).id){
                                node.networkQueues.at(z).q.push(packet);
                            }
                        }
                    }

                }
                temp++;
            }
        }

    }

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
    ///Initilize the Flow and Subflows
    for(int k=0;k<vectorNodes.size();k++){
        vectorNodes.at(k).flow.transDestFlow.first=k;
        vectorNodes.at(k).flow.transDestFlow.second=randomIntNumber(vectorNodes.size());
        while(k==vectorNodes.at(k).flow.transDestFlow.second){
            vectorNodes.at(k).flow.transDestFlow.second=randomIntNumber(vectorNodes.size());
        }
        cout<<"Node "<<vectorNodes.at(k).id<<"has a flow :("<<k<<","<<vectorNodes.at(k).flow.transDestFlow.second<<")"<<"\n";
        unsigned short int numSubFlow=randomIntNumber(maxNumberSubflow);
        if(numSubFlow==0){
            numSubFlow=1;
        }
        vectorNodes.at(k).flow.numberSubflow=numSubFlow;
        cout<<"numSubflow node "<<k<<": is : "<<numSubFlow<<"\n";
        for(int j=0;j<numSubFlow;j++){
            queue<Packet> subFlowQueue;
            vectorNodes.at(k).flow.subflows.push_back(subFlowQueue);
        }


    }
    ///Initilize Routing Table]
    //inistalRoutingTable(nodeNumber);
    ///Show Main Route Table
    //showMainRoutTbale();

    //for(int i=0;i<vectorNodes.size();i++){
      //  initilizeRoutingTable(vectorNodes.at(i),vectorNodes);
    //}
    //showRoutingTable(vectorNodes);

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


                }

        }

    }


    cout<<"-------------------------------------"<<"\n";
    if(numInterface==0){
        cout<<"We have no interface !"<<"\n";
    }



    forwardingTable=distanceVectorAlgorithm(interfaceVector,vectorNodes);
    cout<<"Printing Forwarding Table"<<"\n";
    //printTwoDimentionalVector(forwardingTable);
    cout<<"Graph Bridge is :"<<"\n";
    graph.bridge();


   ///Transferring Packets in each Time Slot
    for(long r=0;r<numberOfTimeSlot;r++){
        for(int t=0;t<vectorNodes.size();t++){
            generatePacket(numPackets,vectorNodes.at(t),r);
            vectorNodes.at(t).flow.flowSize=vectorNodes.at(t).flow.mainFlow.size();
        }

        ///Transferring Packets to the Subflows;
        for(int n=0;n<vectorNodes.size();n++){
                transportToSubflow(vectorNodes.at(n));
        }

        ///Tansport the Packets to the Network Layer Queue
        for(int p=0;p<vectorNodes.size();p++){
            transportToNetwork(vectorNodes.at(p));

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
                            sendPacket(interfaceVector.at(h),vectorNodes.at(recieverIndex),r);
                            intfcTime.at(h)=intfcTime.at(h)-1;
                            }
                        else{
                            int temp=interfaceVector.at(h).interfaceQueue.front().percentSend;
                            sendPacket(interfaceVector.at(h),vectorNodes.at(recieverIndex),r);
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
                            sendPacket(interfaceVector.at(h),vectorNodes.at(recieverIndex),r);
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

    }
    for(int b=0;b<vectorNodes.size();b++){
        cout<<"Recive Queue in Node: "<<b<<" is: "<<vectorNodes.at(b).reciveQueue.size()<<"\n";
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
}
void sendPacket(Interface &interface,Node &node,long arriveTimeSlot){
    if(interface.interfaceQueue.front().destinationAddress==node.id){
        interface.interfaceQueue.front().percentSend=0;
        interface.interfaceQueue.front().arrivaltime=timeSlot*(arriveTimeSlot-interface.interfaceQueue.front().arrivaltime);
        node.reciveQueue.push(interface.interfaceQueue.front());
        interface.interfaceQueue.pop();
    }
    else{
        Packet packet;
        packet=interface.interfaceQueue.front();
        interface.interfaceQueue.pop();
        int source=node.id;
        int destination=packet.destinationAddress;
        int forwarded=forwardingTable.at(source).at(destination);
        for(int i=0;i<node.networkQueues.size();i++){
            if(node.networkQueues.at(i).sourceDestQueue.first==node.id&&
               node.networkQueues.at(i).sourceDestQueue.second==destination){
                node.networkQueues.at(i).q.push(packet);

               }
        }


    }

}
void transportToSubflow(Node &node){

    unsigned long sizeFlow=node.flow.mainFlow.size();
    unsigned long counter=0;
    while(counter<sizeFlow){
        for(int i=0;i<node.flow.numberSubflow;i++){
            if(!node.flow.mainFlow.empty()){
                node.flow.subflows.at(i).push(node.flow.mainFlow.front());
                node.flow.mainFlow.pop();
                counter++;
            }

        }
    }
    for(int j=0;j<node.flow.numberSubflow;j++){
        cout<<"Node "<<node.id<<" subflow "<<j<<" has "<<node.flow.subflows.at(j).size()<<" Packets"<<"\n";
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
/*void updateRoutingTable(Node &transmitter,Node &reciever,int cost,vector<Node> &vectorNodes){
    for(int i=0;i<vectorNodes.size();i++){
            if(transmitter.routingTable.at(i).at(0)==reciever.id){
                transmitter.routingTable.at(i).at(1)=cost;
                transmitter.routingTable.at(i).at(2)=transmitter.id;
            }
    }

}*/
/*void initilizeRoutingTable(Node &node,vector<Node> vectorNodes){
    vector<int> temp;
    for(int j=0;j<vectorNodes.size();j++){
        temp.push_back(j);
    }
    node.routingTable.push_back(temp);
    int valuee=0;
    for(int z=0;z<2;z++){
        vector<int> temp1;
        for(int i=0;i<temp.size();i++){
            temp1.push_back(valuee);
        }
        node.routingTable.push_back(temp1);
        valuee--;
    }
}
void showRoutingTable(vector<Node> vectorNodes){
    for(int i=0;i<vectorNodes.size();i++){
        cout<<"Node "<<i<<"\n";
        for(int j=0;j<vectorNodes.size();j++){
            cout<<vectorNodes.at(i).routingTable.at(j).at(0)<<"  ";
            cout<<vectorNodes.at(i).routingTable.at(j).at(1)<<"  ";
            cout<<vectorNodes.at(i).routingTable.at(j).at(2)<<"\n";

        }
        cout<<"----------------------------------------------"<<"\n";
    }
}
void inistalRoutingTable(unsigned int nodeNumber){
    for(int i=0;i<nodeNumber;i++){
        vector<int> temp;
        for(int j=0;j<nodeNumber;j++){
            temp.push_back(0);
        }
        routingTable.push_back(temp);
    }
}
void showMainRoutTbale(){
    cout<<"Main Routing Table "<<"\n";
    for(int i=0;i<nodeNumber;i++){
        for(int j=0;j<nodeNumber;j++){
            cout<<routingTable.at(i).at(j)<<" ";
        }
        cout<<"\n";
    }
}
void calculateBelmanFordAlgorithm(Node &node,vector<Node> &vectorNodes){
    node.routingTable.at(1).at(node.id)=-1;
    //node.belmanFordWight=0;
    for(int i=0;i<node.strightNeighbors.size();i++){
        int r=node.strightNeighbors.at(i).first;
//        vectorNodes.at(r).belmanFordWight=node.belmanFordWight+1;
//        node.routingTable.at(1).at(r)=vectorNodes.at(r).belmanFordWight;
  //      node.routingTable.at(2).at(r)=node.id;
    }
}
void bellmanFordAlgorithm(int source,vector<Interface> interfaceVector){
        int v=nodeNumber;
        int e=interfaceVector.size();
        int dist[v];
        for(int i=0;i<v;i++){
            dist[i]=INT_MAX;
        }
        dist[source]=0;
        for(int i=1;i<=v-1;i++){
            for(int j=0;j<e;j++){
                int u=interfaceVector.at(j).connectedNode.first;
                int v=interfaceVector.at(j).connectedNode.second;
                int weight=interfaceVector.at(j).weight;
                if(dist[u]!=INT_MAX&&dist[u]+weight<dist[v]){
                    dist[v]=dist[u]+weight;
                }
            }
        }
        ///For detecting the negative Cycle;
        for(int i=0;i<interfaceVector.size();i++){
            int u=interfaceVector.at(i).connectedNode.first;
            int v=interfaceVector.at(i).connectedNode.second;
            int weight=interfaceVector.at(i).weight;
            if(dist[u]!=INT_MAX && dist[u]+weight<dist[v]){
                cout<<"Graph Contain negative Weight Cycle"<<"\n";

            }
        }

}*/
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
                                via[i][j]=k;
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
/*void Graph::connectedComponents()
{
    vector<vector<int>> connections;
    // Mark all the vertices as not visited
    bool *visited = new bool[V];
    for(int v = 0; v < V; v++)
        visited[v] = false;

    for (int v=0; v<V; v++)
    {
        if (visited[v] == false)
        {
            // print all reachable vertices
            // from v
            vector<int> temp;
            connections.push_back(DFSUtil(v, visited));


            cout << "\n";
        }
    }
}
vector<int> Graph::DFSUtil(int v, bool visited[])
{
    // Mark the current node as visited and print it
    visited[v] = true;
    vector<int> dfs;
    cout << v << " ";
    dfs.push_back(v);

    // Recur for all the vertices
    // adjacent to this vertex
    vector<int>::iterator i;
    for(i = adj[v].begin(); i != adj[v].end(); ++i)
        if(!visited[*i])
            DFSUtil(*i, visited);
    return dfs;
}
Graph::Graph(int V)
{
    this->V = V;
    adj = new vector<int>[V];
}

// method to add an undirected edge
void Graph::addEdge(int v, int w)
{
    adj[v].push_back(w);
    adj[w].push_back(v);
}
void printTwoDimentionalVector(vector<vector<int>> table){
    for(int i=0;i<forwardingTable.size();i++){
        for(int j=0;j<forwardingTable.at(i).size();j++){
            cout<<forwardingTable.at(i).at(j)<<" ";
        }
        cout<<"\n";
    }
    cout<<"\n";

}*/



int main()
{

    srand(time(NULL));
    randomTopology();
    return 0;
}
