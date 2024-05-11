#include <iostream>
#include <fstream>
#include <vector>
#include <deque>
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <numeric>
#include <ctime>
#include <random>
#include <sstream>
#include <windows.h>
#include <string>

using namespace std;

template<typename T>
void printElement(const T& element) {
	std::cout << element;
}

void printElement(const std::string& element) {
	std::cout << '\"' << element << '\"';
}

template<typename T>
void printVector(const T& vec, int depth = 0) {
	std::cout << "[";
	bool first = true;
	for (const auto& item : vec) {
		if (!first) {
			std::cout << ", ";
		}
		if constexpr (std::is_class_v<typename T::value_type>) {
			printVector(item, depth + 1);
		}
		else {
			printElement(item);
		}
		first = false;
	}
	std::cout << "]";
}

class RandGen {
public:
	static double normValue(double mean, double stddev, std::mt19937& gen) {
		std::normal_distribution<> dist(mean, stddev);
		return dist(gen);
	}

	static double uniformProb(double lower, double upper, std::mt19937& gen) {
		std::uniform_real_distribution<double> dist(lower, upper);
		return dist(gen);
	}

	static int randInt(int start, int end, std::mt19937& gen) {
		std::uniform_int_distribution<> dist(start, end);
		return dist(gen);
	}
};


// Generates a random Directed Acyclic Graph
vector<vector<int>> GenerateRandomDAG(int numVertices, std::mt19937& gen) {
	int expectedEdges = static_cast<int>(numVertices * log(log(numVertices))); // 期望边数

	// Initialize vertices with values [1, numVertices]
	vector<int> vertices(numVertices);
	for (int i = 0; i < numVertices; ++i) {
		vertices[i] = i + 1;
	}

	// Randomly shuffle vertices to create random edges
	std::shuffle(vertices.begin(), vertices.end(), gen);

	vector<int> neighboringVertices;  // Temporary storage for neighboring vertices
	vector<int> currentAndNeighbors;  // Current vertex and its neighbors
	vector<vector<int>> dag(numVertices);  // The resulting Directed Acyclic Graph

	for (int i = 0; i < numVertices; ++i) {
		auto iter = find(vertices.begin(), vertices.end(), i + 1);

		// Check for neighbor vertices
		if (iter != vertices.begin() && *(iter - 1) < i + 1) {
			neighboringVertices.push_back(*(iter - 1));
		}
		if (iter != vertices.end() - 1 && *(iter + 1) < i + 1) {
			neighboringVertices.push_back(*(iter + 1));
		}

		// Sort and add vertices and its neighbors to the dag
		sort(neighboringVertices.begin(), neighboringVertices.end());
		currentAndNeighbors.push_back(i + 1);
		currentAndNeighbors.insert(currentAndNeighbors.end(), neighboringVertices.begin(), neighboringVertices.end());
		dag[i] = currentAndNeighbors;

		neighboringVertices.clear();
		currentAndNeighbors.clear();
	}

	// Add additional edges to ensure there are expectedEdges edges in the DAG
	int edgeCount = numVertices - 1;
	while (edgeCount < expectedEdges) {
		int randomChild = RandGen::randInt(1, numVertices, gen);
		int randomParent = RandGen::randInt(1, randomChild, gen);

		// Check if the edge already exists
		if (find(dag[randomChild - 1].begin(), dag[randomChild - 1].end(), randomParent) == dag[randomChild - 1].end()) {
			auto iter = dag[randomChild - 1].begin() + 1;

			while (iter != dag[randomChild - 1].end() && *iter <= randomParent) {
				++iter;
			}

			// Insert the edge at the correct position
			dag[randomChild - 1].insert(iter, randomParent);
			++edgeCount;
		}
	}

	return dag;
}


class Bayesnet{
public:
	Bayesnet(){cout<<"the constructor of Bayesnet()"<<endl;};
	~Bayesnet(){cout<<"the destructor of Bayesnet()"<<endl;};
	void print();
	//根据顶点数和连续变量比例以及证据变量个数，对贝叶斯网中的各个元素写值
    void inBayesnet(int Ver, double cpro, int evidnum, std::mt19937& gen);
	void Print();     //输出读入的图
protected:
   vector<vector<int> > verFa;    //V个顶点及其父节点集，每一个是贝叶斯网的基本结构——家族,父集从小号到大号存储
   vector<vector<int> > verChi;   //verChi[k]表示顶点k+1及其子结点的集合
   vector<vector<int> > verdisFa; //V个顶点的家族中离散变量集
   vector<vector<int> > verconFa; //V个顶点的家族中连续变量集
   vector<vector<int> > state;    //V个顶点及其父节点集中具体的离散状态列表,跟家族中离散变量存取顺序有关
   vector<vector<double> > prob;  //V个顶点分布信息列表,跟家族中连续变量的存取顺序有关
   vector<vector<double> > probcon;

   vector<int> distrnum;  //各个顶点拥有的分布个数，由家族中离散变量取几个值决定
   vector<int> dcver;     //描述各个顶点是离散的，还是连续的
   vector<int> valnum;    //描述离散变量的取值个数，连续变量的取值个数记为1

   	vector<int> elabel;           //证据标签，elabel[k]=0表示变量k+1不是证据，elabel[k]=1表示变量k+1是证据
    vector<int> disevid;          //离散证据变量集
	vector<int> disevidsta;       //离散证据变量的状态
	vector<int> conevid;          //连续证据变量集
	vector<double> conevidval;    //连续证据变量的取值


   int V;                //顶点个数
};


void Bayesnet::inBayesnet(int Ver, double cpro, int evidnum, std::mt19937& gen){
       //cout<<"顶点数为 "<<Ver<<" 连续变量比例为 "<<cpro<<" 证据变量个数为 "<<evidnum<<" 的随机贝叶斯网:"<<endl;
       int i,j,m,n,k,l;
       vector<int> vec,ve;
       vector<double> vec1,vec2,vec3;
       V=Ver;
	   verFa = GenerateRandomDAG(V, gen);
       int conv=static_cast<int> (V*cpro);  //连续变量个数
       int disv=V-conv;                     //离散变量个数,让离散变量排在前面

       verChi.clear();
       verdisFa.clear();
       verconFa.clear();
       state.clear();
       prob.clear();
	   probcon.clear();

       distrnum.clear();
       dcver.clear();
       valnum.clear();
       //让标号小的变量是离散变量，标号大的变量是连续变量
       for(i=0;i<V;i++){if(i<disv) dcver.push_back(-1);else dcver.push_back(-2);}
       //这里让离散变量的取值个数都为2
       for(i=0;i<V;i++){if(i<disv) valnum.push_back(2);else valnum.push_back(1);}

       for(k=0;k<V;k++){
	      vec.clear();
		  vec.push_back(k+1);
		  for(l=k+1;l<V;l++){
			  if((find(verFa[l].begin(),verFa[l].end(),k+1))!=verFa[l].end()){
		          vec.push_back(l+1);
			  }
		  }
		  verChi.push_back(vec);
	   }

	   for(k=0;k<V;k++){
		    vec.clear();
		    ve.clear();
		    for(l=0;l<verFa[k].size();l++){
		       if(dcver[verFa[k][l]-1]==-1) vec.push_back(verFa[k][l]);
		       else	   ve.push_back(verFa[k][l]);
		    }
		    verdisFa.push_back(vec);
		    verconFa.push_back(ve);
       }


	   for(k=0;k<V;k++){
		   n=1;
		   for(l=0;l<verFa[k].size();l++){
			    n=n*valnum[verFa[k][l]-1];
		   }
		   distrnum.push_back(n);
	   }

	   for(i=0;i<V;i++){
          ve.clear();
          for(j=0;j<distrnum[i];j++){
              k=j;
              for(l=0;l<verdisFa[i].size();l++){
                  m=k%valnum[verdisFa[i][l]-1];
                  ve.push_back(m+1);
                  k=k/valnum[verdisFa[i][l]-1];
              }
          }
          state.push_back(ve);
       }


       double acc,pr;
	   for (i = 0; i < V; i++) {
		   vec1.clear();
		   vec3.clear();
		   if (dcver[i] == -1) { //离散变量概率的产生机制
			   vec2.clear();
			   for (j = 0; j < distrnum[i]; j++) {
				   vec2.push_back(RandGen::uniformProb(0, 1, gen));
			   }
			   acc = accumulate(vec2.begin(), vec2.end(), 0.0); //计算累积总和

			   // 使用 std::transform 填充 vec1
			   vec1.reserve(vec2.size()); // 预分配所需空间
			   std::transform(vec2.begin(), vec2.end(), std::back_inserter(vec1), [acc](double x) { return x / acc; });

			   vec3 = vec1;

			   for (j = 0; j < vec3.size(); j++) {
				   k = j;
				   if (k % valnum[i] == 0) {
					   acc = 0;
					   for (int jj = 0; jj < valnum[i]; jj++) { acc += vec3[k + jj]; }
					   for (int jj = 0; jj < valnum[i]; jj++) { vec3[k + jj] = vec3[k + jj] / acc; }
				   }
			   }
		   }
		   else {  //连续变量概率的产生机制
			   for (j = 0; j < distrnum[i]; j++) {
				   for (k = 0; k < verconFa[i].size() + 1; k++) {
					   pr = RandGen::uniformProb(1e-10, 1.0, gen);
					   vec1.push_back(pr);
				   }
			   }

			   vec3 = vec1;

		   }
		   prob.push_back(vec1);
		   probcon.push_back(vec3);
       }

    elabel.clear();
   	for(i=0;i<V;i++){elabel.push_back(0);}
    disevid.clear();
	disevidsta.clear();
	conevid.clear();
    conevidval.clear();
    vector<int> tv;
    tv.clear();
    for(i=0;i<V;i++){tv.push_back(i+1); }//设置顶点集
	std::shuffle(tv.begin(), tv.end(), gen); //随机产生证据变量集

    for(j=0;j<evidnum;j++){
       if(dcver[tv[j]-1]==-1){elabel[tv[j]-1]=1;disevid.push_back(tv[j]);}
       else{elabel[tv[j]-1]=1;conevid.push_back(tv[j]);}
    }
    sort(disevid.begin(),disevid.end());
    sort(conevid.begin(),conevid.end());
	for (i = 0; i < disevid.size(); i++) { disevidsta.push_back(RandGen::randInt(1, 2, gen)); }
	for (i = 0; i < conevid.size(); i++) {conevidval.push_back(RandGen::uniformProb(1e-10, 1.0, gen));}
}


void Bayesnet::print(){
     int i,j;
     for(i=0;i<verFa.size();i++){
         for(j=0;j<verFa[i].size();j++){
             //cout<<verFa[i][j]<<" ";
         }//cout<<endl;
     }//cout<<endl;

}


class MoralG : public Bayesnet{     //道义图MoralG,是从Bayesnet继承而来
public:
       MoralG(){cout<<"the constructor of MoralG()"<<endl;};
	   void inMoralG();
	   ~MoralG(){cout<<"the deconstructor of MoralG()"<<endl;};
protected:
       vector<vector<int> > moradj; //moradj[k]表示道义图中顶点k+1及与k+1相邻的顶点集
};

class StarMG : public MoralG{       //星图是从道义图继承而来
public:
       StarMG(){cout<<"the constructor of StarMG()"<<endl;};
       void inStarMG();
	   ~StarMG(){cout<<"the deconstructor of StarMG()"<<endl;};
protected:
	   vector<vector<int> > staradj;//加入星节点V+1,staradj[V]是星节点及其相邻的顶点
};

class OperonG: public StarMG{       //图上的操作,用来构建道义图的极小m三角化图中团的junction tree
public:
	 OperonG(){cout<<"the constructor of OperonG()"<<endl;};
	 void inOperonG();
	 void Snumbering();                 //对星图使用mcs-m算法排序
     int Compare();                     //比较权重大小，返回权重最大的点
	 int Updateweight(int k, int l);    //星图中顶点k被编号的时候，未编号顶点l的权重是否增加

	 void Morsmcs();           //对道义图的极小M三角化图使用Smcs算法
     int Mcompare();           //比较权重大小，返回权重最大的点，如果离散点与连续点权重一样最大，返回离散点
	 void Mupdate(int k);      //顶点k被编号时，更新它的未编号的邻点的权重

     int Complete(vector<int> vec);      //检验点集在道义图的极小M三角化图中是否完全

	 void Constructtree();        //构建道义图的极小m三角化图中团的junction tree
	 int Mverpos();               //记录下点与团的位置的关系,并返回道义图的极小m三角化图中团的个数

	 ~OperonG(){cout<<" the destructor of OperonG: "<<endl;};
protected:
	 vector<int> sweight;           //星图中顶点的当前权重
	 vector<int> stweight;          //星图中顶点的临时权重
	 vector<int> slabel;            //星图中顶点是(1)否(0)已经被编号
	 vector<int> snumbering;        //星图中顶点的编号，snumbering[k]表示顶点k+1的编号
	 vector<int> svisit;            //星图中顶点是(1)否(0)被访问过
	 vector<vector<int> > striadj;  //striadj[k]表示在星图的极小三角化图中点k+1及k+1的相邻点集

     vector<vector<int> > mtriadj;  //mtriadj[k]表示在道义图的极小M三角化图中点k+1及k+1的相邻点集
	 vector<vector<int> > mmadj;    //mmadj[k]表示道义图的极小m三角化图中点k+1及k+1的单调邻集
	 vector<int> mnumbering;        //mnumbering[k]表示点k+1的编号
	 deque<int> morderver;          //按顺序编号的顶点列, morderver[k]表示编号为k+1的点
	 vector<int> mlabel;            //mlabel[k]表示顶点k+1是(1)否(0)已经被编号
	 vector<int> mweight;           //道义图中顶点的当前权重


	 deque<int> msignver;            //msignver[k]表示一个点，这个点标识了排号为k+1的团，即它的单调邻集是分解子
	 deque<int> mcliquever;          //mcliquever[k]表示一个点，这个点及其单调邻集形成排号为k+1的团
     vector<vector<int> > cliques;   //团的列表,cliques[k]表示排号为k+1的团

	 vector<int> mverpos;            //表示使用算法SMCS后在D-numbering下各个点被分配给的那个唯一的团的排号，mverpos[k]表示点k+1所在团的排号
	 vector<int> cliquefa;           //cliquefa[k]表示一个排号，它是团树中排号为k+1的团的父团的排号
	 vector<vector<int> > adjcliques;//adjcliques[k]表示一个排号集，它是团树中与排号为k+1相邻的所有团的排号
	 int cliquenum;                  //团的个数
};

void OperonG::inOperonG(){
   sweight.clear();
   stweight.clear();
   slabel.clear();
   snumbering.clear();
   svisit.clear();
   striadj.clear();

   mtriadj.clear();
   mnumbering.clear();
   mlabel.clear();
   mweight.clear();
   morderver.clear();
   mverpos.clear();

   int i;
   vector<int> vec;
   for(i=0;i<V+1;i++){
      sweight.push_back(0);
	  stweight.push_back(0);
	  slabel.push_back(0);
	  snumbering.push_back(0);
	  svisit.push_back(0);
	  striadj.push_back(staradj[i]);
   }
   for(i=0;i<V;i++){
      mtriadj.push_back(moradj[i]);
	  mnumbering.push_back(0);
	  mlabel.push_back(0);
	  mweight.push_back(0);
	  morderver.push_back(0);
	  mverpos.push_back(0);
   }
   Snumbering();
   Morsmcs();
   Constructtree();
}


class Dismes{
public:
	vector<int> dver;            //离散变量的家族，有些变量可能是证据变量,进入Composterior()后dver仅仅表示一个离散变量集
	vector<vector<int> > dversta;//dversta[k]表示离散变量集dver第k+1个离散状态,是取定dversta中证据变量状态以后的状态空间
	vector<double> dverp;        //dverp[k]表示函数在第k+1个离散状态下的取值,进入Composterior()后dverp仅仅表示一个实数集
};

class Conmes{
public:
	vector<int> cverFa;              //连续变量的家族
    vector<vector<int> > cverFasta;  //连续变量的父集中离散变量的状态列表,是取定离散证据变量状态以后的状态列表
	vector<vector<double> > cverFaf; //对应离散状态列表的条件分布中均值系数及方差列表
};

class Message{
public:
	vector<Conmes> cmes; //团之间传递的连续变量的信息,记录的是条件密度函数信息，每个条件密度函数都是一个信息
    vector<Dismes> dmes; //团之间传递的离散变量的信息，记录的是条件概率函数信息，每个条件概率函数都是一个信息
};


class Dispost{
public:
	int disver;              //顶点
	vector<int> state;    //state[k]表示顶点的第k+1个状态
	vector<double> dispr;    //pr[k]表示顶点在第k+1个状态下的后验概率
};

class Conpost{
public:
    int conver;                     //连续顶点
	vector<double> stap;            //这个点的后验分布是若干个状态下的混合分布，stap[k]表示顶点在第k+1个状态下的概率
    vector<vector<double> > converf;  //conf[k]表示顶点在第k+1个状态下的后验均值和方差。
};


class Propagation: public OperonG{
public:
    Propagation(){cout<<"the constructor of Propagation()"<<endl;};
	void inPropagation();

    void Messagecol();                             //信息收集函数
	Message Cmess(int k);                          //返回值是：收集信息时，团k传向父团的信息potential
	Message Proj(Message & mes, int k, int l);     //将一个团k的信息投影到团l上，返回值也是信息potential,mes表示团k上的当前信息potential
	vector<int> Sep(int k,int l);                  //返回值是：团k和团l的分离子
    Message Projmes(Message & mes, int k, int l); //返回值是：传递的部分信息,在这个函数中使用exchange消去变量
//    Message Elimbarren(vector<int> & RN, Message & mes); //消去非证据变量集RN中的barren变量，并在结构中孤立这些变量,使用这个函数前需要用mes将结构重置
    Message Elimbarren(vector<int> & TAR, vector<int> & RN, Message & mes);

    //返回值是：从点集S出发能够dreach的顶点集(不包含证据变量)，vec是给定的局部证据变量集。
    vector<int> Dreach(vector<int> & S, vector<int> & vec);
    //vec是给定的局部证据变量集，当l是一个碰撞点，是否接受l进一步考虑d连通，返回值是1表示接受，返回值是0表似乎不接受。
    int Acceptl(int l, vector<int> & vec);
    //返回l的后代点中是否与vec中点相交，1表示相交，0表示不相交。
    vector<int> Des(int k);
    //vec是给定的局部证据变量集,由子结点到父结点的边kl被允许记录后，进一步进行判断l和其父或子结点形成的边是否允许被记录。
    void Dreachf(int k, int l, vector<int> & vec);
    //vec是给定的局部证据变量集,由父结点到子结点的边kl被允许记录后，进一步进行判断l和其父子结点形成的边是否允许被记录。
    void Dreachc(int k, int l, vector<int> & vec);
    //返回d分离后用来exchange的信息,dsep是被d分离掉的部分，从mes中去掉所有与desp有关的信息，
    Message Exemess(vector<int> & rels, Message & mes); //即当mes中某个potential涉及的变量集与rels相交非空，则保留这个potential


	int Pos(int k, vector<int> & vec);           //返回值是：变量k在vec中的位置


	vector<int> Remainpos(int k, vector<vector<int> > & vecs);
	//返回变量k的家族中离散变量集在状态vecs下，给定离散证据取值之后，在离散证据变量处一致的状态在vecs中的位置
	vector<vector<int> > Remainsta(int k);
	//返回以变量k为首点的条件密度，在给定离散证据变量取值后，涉及的离散状态
	vector<vector<double> > Remainpr(int k);
	// 返回以变量k为首点的条件密度，在给定离散证据变量取值后，涉及的离散状态对应得密度函数


	vector<vector<int> > Tstate(int k);   //返回点k的动态离散状态列表（考虑父集中的离散证据变量已经给定）
	vector<vector<double> > Tprob(int k); //返回点k的动态条件分布列表（考虑父集中的离散证据变量已经给定）
    int Tpos(int k,int l);  //连续变量k,l,k是l的父点或l自身，返回k在tverconFa[l-1]中的位置
	int Posmess(int k, Message & mes); //若k是连续变量，返回连续变量k在信息mes中连续信息mes.cmes中的位置
	//若k是离散变量，返回离散变量在离散信息mes.dmes中的位置


	vector<Conmes> Exchange(int k, Conmes & conk, int l, Conmes & conl);  //返回exchange后的点k和点l的连续信息
	//在函数及其状态中，对连续变量k及其在tverChi中拓扑序最小的连续变量l进行exchange
    vector<Dismes> Exchange(int k, Dismes & disk, int l, Dismes & disl);  //返回exchange后的点k和点l的离散信息
    //在函数及其状态中，对离散变量k及其在tverChi中拓扑序最小的离散变量l进行exchange
    vector<int> Commonpos(vector<int> & vec, vector<int> & vec_val,
                          vector<int> & ved, vector<vector<int> > & vedsta); //以离散变量集ved为主，
    //返回变量集vec在取值vec_val后，vec与ved的交集上共同的状态在ved的状态中的位置，其中vecsta与vedsta
    //为vec和ved对应的状态。
    vector<int> Commpos(vector<int> & vec, vector<int> & vec_val,  vector<int> & ved);

 	Message Exchange(int k, Message & mes);//在函数上对变量k使用exchange运算,并返回运算后的信息,使用这个函数前需要用mes将结构重置
    vector<vector<int> > EDisstate(vector<int> & vec);//返回离散变量集vec在证据变量取定后的离散状态列表
    void Localstru(Message & mes);   //局部结构化函数，通过信息mes来改写mes涉及的变量的局部结构
    int Domnum(int k);        //返回，为了exchange掉变量k，涉及的离散变量的个数，即k的离散父点，k的子结点的离散父点


	void Messagedis();          //信息分发函数
	Message Dmess(int k);                     //返回值是：分发信息时，团k的父团传向团k的信息potential
    vector<Dismes> Dgenbycv(vector<int> & vec);  //返回值是：由条件密度仅是离散变量函数的连续首点集vec产生的离散信息
    //这个在最后涉及计算后验概率的时候将会有用


    Conpost Comconp(int k, Message & mes);//计算连续点k的后验分布
    Dispost Comdisp(int k, Message & mes);//计算离散点k的后验分布
	void Composterior();        //计算每一个点的后验分布,此时使用VE与VR消元没有区别，因为VE消元可能无法避免，所以，这里对离散变量建议使用VE消元
    double Multip(vector<int> & vec, vector<int> & vecs, vector<Dismes> & dmes); //返回值是：在涉及加和的离散函数信息dmes中，
    //考虑dmes中所有的离散函数带有的离散变量全体vec在固定状态vecs下，对离散函数的取值进行乘积


    double Multiplep(vector<int> inpos, vector<Dismes> & dmes);
    vector<int> Inpos(vector<vector<int> > & allsta, vector<Dismes> & dmes);

    Message Elimdisone(int k, Message & mess); //返回消除离散变量k后留下的信息
    vector<int> Allv(Message & mes);  //返回信息mes中涉及的所有点。
    vector<int> Setminus(vector<int> & allv, vector<int> & parv); //返回allv去掉parv后的点集。

    vector<vector<int> > Disstru(vector<Dismes> & dmes); //返回dmes中非证据离散变量的无向结构
    vector<int> Elimordering(vector<int> & comdisv, vector<Dismes> & dmes);//comdisv中点完全化，返回dmes中的变量集的消元序
    int Maxpoint(vector<int> & comdisv, vector<int> & elimv);//返回elimv中权重最大的点，
    //当有多个权重最大的点且comdisv中有点权重最大时优先返回comdisv中点,comdisv是elimv的子集
    int Changew(int k, int l, vector<int> & elimv);   //考虑是否更新权重


	void Totalmess();   //计算各个团上的总信息，包含原始信息，由父团分发过来的信息，以及由子团传来的信息
    void Output(vector<Dispost> & dis, vector<Conpost> & con);  //以文本形式输出最终的结果
    void Outputstru(int Ver, double cpro, int evidnum, int net_id, std::mt19937& gen);
	void OutputNet(int Ver, double cpro, int evidnum, int net_id, std::mt19937& gen);
	void OutputData(int Ver, double cpro, int evidnum, int n_sample, int net_id, std::mt19937& gen);
	~Propagation();
private:
	vector<Message> omess;   //omess[k]表示团k+1上的原始信息potential
	vector<Message> smess1;  //smess1[k]表示团k+1上传向它的父团的信息potential,smess1[0]表示团1上的原始信息potential
	vector<Message> smess2;  //smess2[k]表示由团k+1的父团传向团k+1的信息potential,smess2[0]表示团1上的原始信息potential
    vector<int> meslabel;    //meslabel[k]=0表示团k+1还没有从它的父团接受信息，meslabel[k]=1表示团k+1已经从它的父团接受信息
    vector<Message> totalmess; //totalmess[k]表示团k+1上所有信息，包含原始信息，从父团传来的信息以及从子团传来的信息

    vector<int> verfapos;            //verfapos[k]表示一个排号，它是包含顶点k+1及其父集的一个团的排号
	vector<vector<int> > cliquefun;  //cliquefun[k]表示一个顶点集，它是排号为k+1的团中被赋于的条件概率,条件密度的首点集

    //以下7个变量集会在exchange中有所改变
	vector<vector<int> > tverChi;     //临时的结点及其子结点集合
	vector<vector<int> > tverFa;      //临时的家族集合
    vector<vector<int> > tverdisFa;   //临时的V个顶点的家族中离散变量集
    vector<vector<int> > tverconFa;   //临时的V个顶点的家族中连续变量集
	vector<vector<int> > tstate;   //临时的离散状态列表
	vector<vector<double> > tprob; //临时的条件分布
	vector<int> tdistrnum;         //临时的各个变量涉及的分布个数

	vector<vector<int> > estate;   //考虑离散证据变量取值后的离散状态列表
	vector<vector<double> > eprob; //考虑离散证据变量取值后的条件分布
	vector<int> edistrnum;         //考虑离散证据变量取值后的各个变量涉及的分布个数，这3个表在Propagation()构建时确定下来

    vector<int> desvisit; //在深度优先遍历后代结点时，记录是否已经遍历了这个顶点，1表示已经遍历，0表示没有遍历
    vector<int> dreach; //记录点是否已经dreach,1表示dreach,0表示没有
    vector<vector<int> > Favisit;//子父边是否访问过，1表示访问过，0表示没有
    vector<vector<int> > Chivisit;//父子边是否访问过，1表示访问过，0表示没有

    vector<int> evid;        //证据变量集
    vector<int> poster;      //非证据变量集，需要计算其中每个点的后验概率

     vector<int> weight;           //顶点的当前权重
	 vector<int> tweight;          //顶点的临时权重
	 vector<int> label;            //顶点是(1)否(0)已经被编号
	 vector<int> visit;            //顶点是(1)否(0)被访问过
	 vector<vector<int> > stru;    //顶点间结构

};


vector<int> Propagation::Dreach(vector<int> & S, vector<int> & vec){
//       cout<<"调用Dreach: "<<endl;
       int i,j,k;
       vector<int> dreachve;
       for(i=0;i<V;i++){dreach[i]=0;}
       Favisit=tverFa; //使用Dreach时，先将Favisit和Chivisit重置
       Chivisit=tverChi;
       for(i=0;i<Favisit.size();i++){for(j=0;j<Favisit[i].size();j++){Favisit[i][j]=0;}}
       for(i=0;i<Chivisit.size();i++){for(j=0;j<Chivisit[i].size();j++){Chivisit[i][j]=0;}}
       for(j=0;j<S.size();j++){
          k=S[j];
          dreach[k-1]=1;//认为S中的点与自身d可达是一种平凡的情况
          for(i=1;i<tverFa[k-1].size();i++){
              //cout<<"dreachf( "<<k<<" , "<<tverFa[k-1][i]<<" )"<<endl;
              Favisit[k-1][i]=1;Dreachf(k,tverFa[k-1][i],vec);
              if(elabel[tverFa[k-1][i]-1]==0) dreach[tverFa[k-1][i]-1]=1;
          }
          for(i=1;i<tverChi[k-1].size();i++){
             // cout<<"dreachc( "<<k<<" , "<<tverChi[k-1][i]<<" )"<<endl;
              Chivisit[k-1][i]=1;Dreachc(k,tverChi[k-1][i],vec);
              if(elabel[tverChi[k-1][i]-1]==0) dreach[tverChi[k-1][i]-1]=1;
          }
       }
       for(i=0;i<V;i++){if(dreach[i]==1) dreachve.push_back(i+1);}
       //cout<<"dreach的顶点"<<endl;
      // for(i=0;i<dreachve.size();i++){cout<<dreachve[i]<<" ";}cout<<endl;
       return dreachve;
}

void Propagation::Dreachf(int k, int l, vector<int> & vec){
        int i;
        if(elabel[l-1]==0){
             for(i=1;i<tverFa[l-1].size();i++){
//             cout<<tverFa[l-1][i]<<" "<<Favisit[l-1][i]<<endl;
               if(Favisit[l-1][i]==0){
//                  cout<<"11111Dreachf( "<<l<<" , "<<tverFa[l-1][i]<<" )"<<endl;
                  Favisit[l-1][i]=1;Dreachf(l,tverFa[l-1][i],vec);
                  if(elabel[tverFa[l-1][i]-1]==0) dreach[tverFa[l-1][i]-1]=1;
               }
             }
             for(i=1;i<tverChi[l-1].size();i++){
//               cout<<tverChi[l-1][i]<<" "<<Chivisit[l-1][i]<<endl;
               if((Chivisit[l-1][i]==0)&&(tverChi[l-1][i]!=k)){
//                  cout<<"1111Dreachc( "<<l<<" , "<<tverChi[l-1][i]<<" )"<<endl;
                  Chivisit[l-1][i]=1;
                  Dreachc(l,tverChi[l-1][i],vec);
                  if(elabel[tverChi[l-1][i]-1]==0) dreach[tverChi[l-1][i]-1]=1;
               }
             }
        }
}

void Propagation::Dreachc(int k, int l, vector<int> & vec){
        int i;
        for(i=1;i<tverFa[l-1].size();i++){
//            cout<<tverFa[l-1][i]<<" "<<Acceptl(l,vec)<<" "<<Favisit[l-1][i]<<endl;
            if((Acceptl(l,vec)==1)&&(Favisit[l-1][i]==0)&&(tverFa[l-1][i]!=k)){
//               cout<<"2222Dreachf( "<<l<<" , "<<tverFa[l-1][i]<<" )"<<endl;
                Favisit[l-1][i]=1;
                Dreachf(l,tverFa[l-1][i],vec);
                if(elabel[tverFa[l-1][i]-1]==0) dreach[tverFa[l-1][i]-1]=1;
            }
        }
        for(i=1;i<tverChi[l-1].size();i++){
//             cout<<tverChi[l-1][i]<<" "<<elabel[l-1]<<" "<<Chivisit[l-1][i]<<endl;
            if((elabel[l-1]==0)&&(Chivisit[l-1][i]==0)){
//                cout<<"2222Dreachc( "<<l<<" , "<<tverChi[l-1][i]<<" )"<<endl;
                Chivisit[l-1][i]=1;Dreachc(l,tverChi[l-1][i],vec);
                if(elabel[tverChi[l-1][i]-1]==0) dreach[tverChi[l-1][i]-1]=1;
            }
        }
}

vector<int> Propagation::Des(int k){
    desvisit[k-1]=1;
    vector<int> vec,vec1;
    vec.push_back(k);
    int i,j;
    for(i=1;i<tverChi[k-1].size();i++){
      if(desvisit[tverChi[k-1][i]-1]==0){
           vec1=Des(tverChi[k-1][i]);
           for(j=0;j<vec1.size();j++){vec.push_back(vec1[j]);}
       }
    }
    return vec;
}

int Propagation::Acceptl(int l, vector<int> & vec){
    int i;
    vector<int> vect; //vect描述l和l的后代集的并集
    if(vec.empty()) return 0;
    else{
        desvisit.clear();
        for(i=0;i<V;i++){desvisit.push_back(0);}
        vect=Des(l);
        for(i=0;i<vect.size();i++){
            if(find(vec.begin(),vec.end(),vect[i])!=vec.end()) return 1;
        }
    }
    return 0;
}




void Propagation::inPropagation(){
	int i,j,q,s,k;
	vector<int> vec1,vec2,ve,vepos,pos1;
	vector<int>::iterator iter;
    tverFa=verFa;
	tverChi=verChi;
	tverdisFa=verdisFa;
	tverconFa=verconFa;

	omess.clear();
    meslabel.clear();

	desvisit.clear();
	dreach.clear();
    verfapos.clear();
    cliquefun.clear();
    evid.clear();
    poster.clear();


    for(i=0;i<V;i++){desvisit.push_back(0);}
    for(i=0;i<V;i++){dreach.push_back(0);}
    Favisit=tverFa;
    Chivisit=tverChi;

	estate=state;
    eprob=prob;
	edistrnum=distrnum;

	tstate=state;
    tprob=prob;
	tdistrnum=distrnum;


    for(i=0;i<V;i++){
		vec1=verFa[i];
		vec2.clear();
		for(j=0;j<vec1.size();j++){vec2.push_back(mnumbering[vec1[j]-1]);}
		iter=min_element(vec2.begin(),vec2.end());       //取出家族中最小的编号
		k=morderver[*iter-1];                            //最小编号对应的顶点
	    verfapos.push_back(mverpos[k-1]);
	}

    //顶点的分配方式需要重写
	for(i=0;i<cliquenum;i++){
		vec1.clear();
		for(j=0;j<verfapos.size();j++){if(verfapos[j]==i+1) vec1.push_back(j+1);}
		cliquefun.push_back(vec1);
	}
	for(i=0;i<V;i++){if(elabel[i]==0) poster.push_back(i+1);else evid.push_back(i+1);}

	Message mes;
	Conmes conm;
	Dismes dism;
	vector<vector<int> > vec,vecs;
	vector<vector<double> > vecp;
	for(i=0;i<cliquenum;i++){
		mes.cmes.clear();
		mes.dmes.clear();
		vec1.clear();
		vec2.clear();
	    for(j=0;j<cliquefun[i].size();j++){if(dcver[cliquefun[i][j]-1]==-2) vec1.push_back(cliquefun[i][j]);else vec2.push_back(cliquefun[i][j]);}
		for(j=0;j<vec1.size();j++){
			conm.cverFaf.clear();
			conm.cverFasta.clear();
			conm.cverFa=verFa[vec1[j]-1];
			vecs=Remainsta(vec1[j]);
			vecp=Remainpr(vec1[j]);
			estate[vec1[j]-1].clear();
			eprob[vec1[j]-1].clear();
			for(q=0;q<vecs.size();q++){for(s=0;s<vecs[q].size();s++){estate[vec1[j]-1].push_back(vecs[q][s]);}}
			for(q=0;q<vecp.size();q++){for(s=0;s<vecp[q].size();s++){eprob[vec1[j]-1].push_back(vecp[q][s]);}}
			edistrnum[vec1[j]-1]=vecs.size();
			for(q=0;q<vecs.size();q++){conm.cverFasta.push_back(vecs[q]);}
			for(q=0;q<vecp.size();q++){conm.cverFaf.push_back(vecp[q]);}
		    mes.cmes.push_back(conm);

		}
		for(j=0;j<vec2.size();j++){
			dism.dverp.clear();
			dism.dversta.clear();
			dism.dver=tverFa[vec2[j]-1];
			vecs=Remainsta(vec2[j]);
			vecp=Remainpr(vec2[j]);
            estate[vec2[j]-1].clear();
			eprob[vec2[j]-1].clear();
			for(q=0;q<vecs.size();q++){for(s=0;s<vecs[q].size();s++){estate[vec2[j]-1].push_back(vecs[q][s]);}}
			for(q=0;q<vecp.size();q++){for(s=0;s<vecp[q].size();s++){eprob[vec2[j]-1].push_back(vecp[q][s]);}}
			edistrnum[vec2[j]-1]=vecs.size();
			for(q=0;q<vecs.size();q++){dism.dversta.push_back(vecs[q]);}
			for(q=0;q<vecp.size();q++){dism.dverp.push_back(vecp[q][0]);}
			mes.dmes.push_back(dism);
		}

		omess.push_back(mes);
	}

	tstate=estate;
    tprob=eprob;
	tdistrnum=edistrnum;
	smess1=omess;
	smess2=omess;
	totalmess=omess;

	for(i=0;i<cliquenum;i++){meslabel.push_back(0);}
}

void Propagation::Messagecol(){
   //cout<<"开始信息收集： "<<endl;
   int i;
   for(i=0;i<adjcliques[0].size();i++){
	   if(adjcliques[0][i]!=cliquefa[0]){
		   //cout<<"团1调用团 "<<adjcliques[0][i]<<" 传向团1的信息"<<endl;
		   Cmess(adjcliques[0][i]);
	   }
   }
   //cout<<"信息收集完毕 "<<endl;
}


void Propagation::Localstru(Message & mes){
	//cout<<"调用Localstru() "<<endl;
	int i,j,m;
    vector<int> vec1;
	vector<vector<int> > vec;         //vec是信息中家族的集合
	for(i=0;i<mes.cmes.size();i++){
		vec.push_back(mes.cmes[i].cverFa);
	}
	for(i=0;i<mes.dmes.size();i++){
        vec.push_back(mes.dmes[i].dver);
    }

	for(i=0;i<V;i++){tverFa[i].clear();tverFa[i].push_back(i+1);}
	for(i=0;i<vec.size();i++){ //确定家族集
	    tverFa[vec[i][0]-1]=vec[i];
	}
	for(i=0;i<V;i++){
		tverdisFa[i].clear();
		tverconFa[i].clear();
		for(j=0;j<tverFa[i].size();j++){
		    if(dcver[tverFa[i][j]-1]==-1) tverdisFa[i].push_back(tverFa[i][j]);
			else tverconFa[i].push_back(tverFa[i][j]);
		}
	}

 	for(i=0;i<V;i++){
	    vec1.clear();
		vec1.push_back(i+1);
		for(j=i+1;j<V;j++){
			if((find(tverFa[j].begin(),tverFa[j].end(),i+1))!=tverFa[j].end()){
		          vec1.push_back(j+1);
			}
		}
		tverChi[i]=vec1;
	}

	vector<Conmes> conmes=mes.cmes;  //利用conmes对连续信息进行重置
	vector<Dismes> dismes=mes.dmes;   //利用dismes对离散信息进行重置
	int conv,disv;
	for(i=0;i<conmes.size();i++){
        conv=conmes[i].cverFa[0];
		tstate[conv-1].clear();
		tprob[conv-1].clear();
		tdistrnum[conv-1]=conmes[i].cverFasta.size();
		for(j=0;j<conmes[i].cverFasta.size();j++){
			for(m=0;m<conmes[i].cverFasta[j].size();m++){tstate[conv-1].push_back(conmes[i].cverFasta[j][m]);}
		}
		for(j=0;j<conmes[i].cverFaf.size();j++){
			for(m=0;m<conmes[i].cverFaf[j].size();m++){tprob[conv-1].push_back(conmes[i].cverFaf[j][m]);}
		}
	}
	for(i=0;i<dismes.size();i++){
        disv=dismes[i].dver[0];
        tstate[disv-1].clear();
        tprob[disv-1].clear();
        tdistrnum[disv-1]=dismes[i].dversta.size();
        for(j=0;j<dismes[i].dversta.size();j++){
              for(m=0;m<dismes[i].dversta[j].size();m++){tstate[disv-1].push_back(dismes[i].dversta[j][m]);}
        }
        for(j=0;j<dismes[i].dverp.size();j++){tprob[disv-1].push_back(dismes[i].dverp[j]);}
    }
	//cout<<"调用Localstru()完毕 "<<endl<<endl;
}


Message Propagation::Exemess(vector<int> & rels, Message & mes){
        int i,j;
        vector<int>::iterator it;
        Message mess;
        for(i=0;i<mes.cmes.size();i++){
            for(j=0;j<mes.cmes[i].cverFa.size();j++){
                it=find(rels.begin(),rels.end(),mes.cmes[i].cverFa[j]);
                if(it!=rels.end()){mess.cmes.push_back(mes.cmes[i]);break;}
            }

        }
        for(i=0;i<mes.dmes.size();i++){
            for(j=0;j<mes.dmes[i].dver.size();j++){
                it=find(rels.begin(),rels.end(),mes.dmes[i].dver[j]);
                if(it!=rels.end()){mess.dmes.push_back(mes.dmes[i]);break;}
            }
        }
        return mess;
}


Message Propagation::Elimbarren(vector<int> & TAR, vector<int> & RN, Message & mes){
     int i,j;
     vector<int> des,barren;
     Message mess;
     for(i=0;i<RN.size();i++){
         desvisit.clear();
         for(j=0;j<V;j++){desvisit.push_back(0);}
         des=Des(RN[i]);
         for(j=0;j<des.size();j++){
             if(elabel[des[j]-1]==1) break;
             if(find(TAR.begin(),TAR.end(),des[j])!=TAR.end()) break;
         }
         if(j==des.size()) barren.push_back(RN[i]);
     }


     for(i=0;i<mes.cmes.size();i++){
         if(find(barren.begin(),barren.end(),mes.cmes[i].cverFa[0])==barren.end())
            mess.cmes.push_back(mes.cmes[i]);
     }
     for(i=0;i<mes.dmes.size();i++){
         if(find(barren.begin(),barren.end(),mes.dmes[i].dver[0])==barren.end())
            mess.dmes.push_back(mes.dmes[i]);
     }
     return mess;
}


Message Propagation::Cmess(int k){ //返回值是：收集信息时，团k传向父团的信息potential
	int i,j;
	Message mes=omess[k-1];       //mes表示团k上的当前信息
	vector<int> vec=adjcliques[k-1];
	Message tmes,cmes;
	if(vec.size()==1){
	   cmes=Projmes(mes,k,cliquefa[k-1]);
	}
    else{
		for(i=0;i<vec.size();i++){
		//	cout<<"团 "<<k<<" 的邻团 "<<vec[i]<<" "<<endl;
			if(vec[i]!=cliquefa[k-1]){
		//		cout<<"团 "<<k<<" 调用团"<<vec[i]<<" 传向团 "<<k<<" 的信息 "<<endl;
			    tmes=Cmess(vec[i]);
			    for(j=0;j<tmes.cmes.size();j++){
			      mes.cmes.push_back(tmes.cmes[j]);
				}
			    for(j=0;j<tmes.dmes.size();j++){
			      mes.dmes.push_back(tmes.dmes[j]);
				}
			}
		}
		cmes=Projmes(mes,k,cliquefa[k-1]);
	}
	smess1[k-1]=cmes;         //记录下信息收集时，团k传向它的父团的信息
    return cmes;
}

Message Propagation::Dmess(int k){        //k应该不为1,即不考虑团1的父团0传向团1的信息
	//cout<<"调用函数Dmess()"<<endl;
	int i,j,m;
	m=cliquefa[k-1];             //团k的父团是m
	Message mes=omess[m-1];      //团m上带有的原始信息

	vector<int> vec=adjcliques[m-1];  //团m的相邻团
	Message tmes,cmes;
//	cout<<"信息分发过程开始："<<endl;
//	cout<<"团 "<<m<<" 向团 "<<k<<" 分发信息"<<endl;
	for(i=0;i<vec.size();i++){
		if((vec[i]!=k)&&(vec[i]!=0)){	//若团m的相邻团vec[i]既不是团0（即假想的团1的父团）又不是团k
//			cout<<"信息分发时"<<"团 "<<m<<" 调用团 "<<vec[i]<<" 传来的信息"<<endl;
			if(cliquefa[vec[i]-1]==m){
//				cout<<"信息分发时"<<"团 "<<m<<" 调用团 "<<vec[i]<<" 传来的信息已经在信息收集时被记录"<<endl;
			    tmes=smess1[vec[i]-1];

			}
			else {
//				cout<<"因为团 "<<vec[i]<<" 是团 "<<m<<" 的父团,进一步调用信息分发"<<endl;
				if(meslabel[m-1]==0){
//					cout<<"团 "<<m<<" 未接受父团传来的信息"<<endl;
					tmes=Dmess(m);
				}
				else{
//					cout<<"团 "<<m<<" 在之前信息分发过程中，已经接受了父团 "<<vec[i]<<" 传来的信息"<<endl;
					tmes=smess2[m-1];
				}

     /*           cout<<"连续信息2: "<<endl;
				for(j=0;j<tmes.cmes.size();j++){
					for(int q=0;q<tmes.cmes[j].cverFa.size();q++){
					    cout<<tmes.cmes[j].cverFa[q]<<" ";
					}cout<<endl;
				}cout<<endl;
				cout<<"离散信息2: "<<endl;
				for(j=0;j<tmes.dmes.size();j++){
					for(int q=0;q<tmes.dmes[j].dver.size();q++){
					    cout<<tmes.dmes[j].dver[q]<<" ";
					}cout<<endl;
				}cout<<endl;*/

			}
		   for(j=0;j<tmes.cmes.size();j++){mes.cmes.push_back(tmes.cmes[j]);}
		   for(j=0;j<tmes.dmes.size();j++){mes.dmes.push_back(tmes.dmes[j]);}
		}
	}

	cmes=Projmes(mes,m,k);
    smess2[k-1]=cmes;    //记录分发信息时，团k从它的父团接受的信息
	meslabel[k-1]=1;
	//cout<<"调用函数Dmess( "<<k<<" )完毕"<<endl;
    return cmes;
}


Message Propagation::Projmes(Message & mes, int k, int l){
	//cout<<"调用函数Projmes()："<<endl;
	//用A表示团k上的当前信息点集，一部分是条件密度的首点，一部分是条件概率的首点
	vector<int> A,R,RE,RN,RNC,RND,RNDE;
	vector<int>::iterator it;
    int i,j,m,num;
	vector<int> S=Sep(k, l);
	A=Allv(mes);R=A;
	for(i=0;i<S.size();i++){it=find(R.begin(),R.end(),S[i]);if(it!=R.end()){R.erase(it);}}
	for(i=0;i<R.size();i++){if(elabel[R[i]-1]==0) RN.push_back(R[i]);else RE.push_back(R[i]);}
	Localstru(mes);
	Message remess=Elimbarren(S,RN,mes),mess,mess1;//remess表示消去barren变量后余下的信息
	Localstru(remess);
	vector<int> relsvec=Dreach(S,RE);//S能够d可达的点集
    mess=Exemess(relsvec,remess);
    Localstru(mess);


    RN.clear();A=Allv(mess);R=A;
   	for(i=0;i<S.size();i++){it=find(R.begin(),R.end(),S[i]);if(it!=R.end()){R.erase(it);}}
	for(i=0;i<R.size();i++){if(elabel[R[i]-1]==0) RN.push_back(R[i]);}
    for(i=0;i<RN.size();i++){if(dcver[RN[i]-1]==-1) RND.push_back(RN[i]);else RNC.push_back(RN[i]);}

    for(i=RNC.size()-1;i>=0;i--){num=Domnum(RNC[i]);if(num<14) mess=Exchange(RNC[i],mess);}
	//RNC已经被从小到大排序，从序号大的点开始exchange消去

	//残差中的非证据离散变量集中没有连续变量作为子结点的离散变量可以被exchange掉
	for(i=0;i<RND.size();i++){
            m=0;
            for(j=1;j<tverChi[RND[i]-1].size();j++){
                      if(dcver[tverChi[RND[i]-1][j]-1]==-2) {m++; break;}
            }
            if(m==0) RNDE.push_back(RND[i]);
    }
	for(i=RNDE.size()-1;i>=0;i--){num=Domnum(RNDE[i]);if(num<14) mess=Exchange(RNDE[i],mess);}
    //cout<<"调用函数Projmes()完毕"<<endl;
    return mess;
}

int Propagation::Domnum(int k){
	int i,j,m;
	vector<int> dom;
	int num;
	for(i=0;i<tverdisFa[k-1].size();i++){if(elabel[tverdisFa[k-1][i]-1]==0)  dom.push_back(tverdisFa[k-1][i]);}
    for(i=0;i<tverChi[k-1].size();i++){
		m=tverChi[k-1][i];
		for(j=0;j<tverdisFa[m-1].size();j++){if(elabel[tverdisFa[m-1][j]-1]==0) dom.push_back(tverdisFa[m-1][j]);}
	}
	sort(dom.begin(),dom.end());
	vector<int>::iterator it=unique(dom.begin(),dom.end());
	dom.erase(it,dom.end());
    num=dom.size();
	return num;
}


vector<vector<int> > Propagation::EDisstate(vector<int> & vec){
         int i,j;
		 int n=1;
		 int m,k,s,t;
		 vector<vector<int> > sta,evista;
		 vector<int> ve,evec,epos,nevec,ne;
         for(i=0;i<vec.size();i++){
			if(elabel[vec[i]-1]==1){evec.push_back(vec[i]);}
			else{nevec.push_back(vec[i]);}
		 }
         for(i=0;i<evec.size();i++){epos.push_back(Pos(evec[i],disevid));}
         for(i=0;i<nevec.size();i++){n=n*valnum[nevec[i]-1];}
         for(i=0;i<n;i++){
			 ve.clear();
			 m=i;
			 for(j=0;j<nevec.size();j++){
				 k=m%valnum[nevec[j]-1];
				 ve.push_back(k+1);
				 m=m/valnum[nevec[j]-1];
			 }
			 sta.push_back(ve);
		 }

		 for(i=0;i<sta.size();i++){
		     ne.clear();
             s=0;
             t=0;
		     for(j=0;j<vec.size();j++){
		         if(elabel[vec[j]-1]==1){ne.push_back(disevidsta[epos[s]-1]);s++;}
		         else{ne.push_back(sta[i][t]);t++;}
		     }
             evista.push_back(ne);
		 }
		 return evista;
}


vector<int> Propagation::Inpos(vector<vector<int> > & allsta, vector<Dismes> & dmes){
    int i,j;
    int pos,power;
	vector<int> inpos;         //利用十进制和多进制的转换关系
	for(i=0;i<allsta.size();i++){
	    pos=0;power=1;
	    for(j=0;j<allsta[i].size();j++){
	        if(elabel[dmes[i].dver[j]-1]==0){
	            pos=pos+(allsta[i][j]-1)*power;
	            power=power*valnum[dmes[i].dver[j]-1];
	        }
	    }
	    inpos.push_back(pos);
	}
    return inpos;
}

vector<int> Propagation::Setminus(vector<int> & allv, vector<int> & parv){
      int i;
      vector<int> remv=allv;
      vector<int>::iterator it;
      for(i=0;i<parv.size();i++){
          it=find(remv.begin(),remv.end(),parv[i]);
          if(it!=remv.end()) remv.erase(it);
      }
      return remv;
}

vector<int> Propagation::Commpos(vector<int> & vec, vector<int> & vec_val,
               vector<int> & ved){
        int i,j;
        int pos,power;
        vector<int> commpos,veds;
        vector<int> remk=Setminus(ved,vec); //remk中最多只有一个元素k
        if(!remk.empty()){
           for(i=0;i<valnum[remk[0]-1];i++){
               veds.clear();
               for(j=0;j<ved.size();j++){
                   if(ved[j]==remk[0]) veds.push_back(i+1);
                   else veds.push_back(vec_val[Pos(ved[j],vec)-1]);
               }
               pos=0;power=1;
               for(j=0;j<ved.size();j++){
                   if(elabel[ved[j]-1]==0){
                      pos=pos+(veds[j]-1)*power;
                      power=power*valnum[ved[j]-1];
                   }
               }
               commpos.push_back(pos+1);
           }
        }
        else{
            veds.clear();
            for(j=0;j<ved.size();j++){veds.push_back(vec_val[Pos(ved[j],vec)-1]);}
            pos=0;power=1;
            for(j=0;j<ved.size();j++){
                if(elabel[ved[j]-1]==0){
                   pos=pos+(veds[j]-1)*power;
                   power=power*valnum[ved[j]-1];
                }

            }
            commpos.push_back(pos+1);
        }
        return commpos;
}


vector<int> Propagation::Commonpos(vector<int> & vec, vector<int> & vec_val,
               vector<int> & ved, vector<vector<int> > & vedsta){
               int i,j,m;
               vector<int> comm1pos,comm2pos,vedstapos,ve;
               for(i=0;i<vec.size();i++){
                   if(find(ved.begin(),ved.end(),vec[i])!=ved.end()){
                       comm1pos.push_back(i+1);
                       comm2pos.push_back(Pos(vec[i],ved));
                       }
               }
               for(i=0;i<vedsta.size();i++){
                    m=0;
                    ve=vedsta[i];
                    for(j=0;j<comm1pos.size();j++){
                          if(vec_val[comm1pos[j]-1]!=ve[comm2pos[j]-1]) {m++;break;}
                    }
                    if(m==0) vedstapos.push_back(i+1);
               }
               return vedstapos;
}




vector<Conmes> Propagation::Exchange(int k, Conmes & conk, int l, Conmes & conl){ //在函数及其状态中，对连续变量k及其在tverChi中拓扑序最小的连续变量l进行exchange
    //cout<<"调用函数Exchange("<<k<<" , "<<l<<")连续变量"<<endl;
	int i,j;
	Conmes Conm;
	vector<Conmes>  Vecconm;
	vector<int> veA,veB;
	vector<double> fAk,fBl,fk,fl;
	vector<int> vecd1,vecd2,vecc1,vecc2;
	vector<vector<int> > st1,st2,st3;     //exchange后k和l有共同的离散状态st3
	vector<vector<double> > pr1,pr2,prk,prl;
	vector<int>::iterator it;
	st1=conk.cverFasta;
	st2=conl.cverFasta;
	pr1=conk.cverFaf;
	pr2=conl.cverFaf;
	vector<int> D(tverdisFa[k-1].size()+tverdisFa[l-1].size());
	vecd1=tverdisFa[k-1];      //暂时固定住k与l的离散父集
	vecd2=tverdisFa[l-1];
	merge(vecd1.begin(),vecd1.end(),vecd2.begin(),vecd2.end(),D.begin()); //将k与l的离散家族集并起来
    sort(D.begin(),D.end());
	it=unique(D.begin(),D.end());
	D.erase(it,D.end());                 //去掉D中重复的点
	st3=EDisstate(D);
    vector<int> C(tverconFa[k-1].size()+tverconFa[l-1].size());
    vecc1=tverconFa[k-1];       //暂时固定住k与l的连续父集
    vecc2=tverconFa[l-1];
	sort(vecc1.begin(),vecc1.end());
	sort(vecc2.begin(),vecc2.end());
	merge(vecc1.begin(),vecc1.end(),vecc2.begin(),vecc2.end(),C.begin()); //将k与l的连续家族集并起来
	sort(C.begin(),C.end());
    it=unique(C.begin(),C.end());
	C.erase(it,C.end());                //去掉C中重复的点
	it=find(C.begin(),C.end(),k);
	if(it!=C.end()) C.erase(it);        //C实际上是exchange后k的连续父结点集
	vector<int> CC=C;//CC是exchange后k与l共同的连续父结点集，即从C中去掉点l可得CC
	it=find(CC.begin(),CC.end(),l);
	if(it!=CC.end()) CC.erase(it);
    int pos=Pos(k,tverconFa[l-1]);  //记录下k在l的连续家族集中的位置
	double cof;
    double a,b;
	for(i=0;i<st3.size();i++){
		veA.clear();
        veB.clear();
		fl.clear();
		fk.clear();
		for(j=0;j<st3[i].size();j++){
			if(find(vecd1.begin(),vecd1.end(),D[j])!=vecd1.end()) veA.push_back(st3[i][j]);
			if(find(vecd2.begin(),vecd2.end(),D[j])!=vecd2.end()) veB.push_back(st3[i][j]);
		}
		for(j=0;j<st1.size();j++){if(equal(veA.begin(),veA.end(),st1[j].begin())){fAk=pr1[j];break;}}
		for(j=0;j<st2.size();j++){if(equal(veB.begin(),veB.end(),st2[j].begin())){fBl=pr2[j];break;}}
		cof=fBl[pos-1];
        fl.push_back(fBl[0]+cof*fAk[0]); //exchange后l的回归常数
		//exchange后k的回归常数
		fk.push_back((fAk[0]*fBl[fBl.size()-1]-fBl[0]*cof*fAk[fAk.size()-1])/(fBl[fBl.size()-1]+cof*cof*fAk[fAk.size()-1]));
		for(j=0;j<CC.size();j++){
			if(find(vecc1.begin(),vecc1.end(),CC[j])!=vecc1.end()) a=fAk[Pos(CC[j],tverconFa[k-1])-1];//a=fAk[Tpos(CC[j],k)-1];
			else a=0;//
			if(find(vecc2.begin(),vecc2.end(),CC[j])!=vecc2.end()) b=fBl[Pos(CC[j],tverconFa[l-1])-1];//b=fBl[Tpos(CC[j],l)-1];
			else b=0; //
            fl.push_back(b+cof*a); //exchange后l关于CC=C\{l}回归系数
			//exchange后k关于CC=C\{l}回归系数
			fk.push_back((a*fBl[fBl.size()-1]-b*cof*fAk[fAk.size()-1])/(fBl[fBl.size()-1]+cof*cof*fAk[fAk.size()-1]));
		}

		fl.push_back(fBl[fBl.size()-1]+cof*cof*fAk[fAk.size()-1]); //exchange后l的方差
		prl.push_back(fl);
		fk.push_back(cof*fAk[fAk.size()-1]/(fBl[fBl.size()-1]+cof*cof*fAk[fAk.size()-1])); //exchange后k关于l的回归系数
		fk.push_back(fAk[fAk.size()-1]*fBl[fBl.size()-1]/(fBl[fBl.size()-1]+cof*cof*fAk[fAk.size()-1])); //exchange后k的方差
		prk.push_back(fk);
	}

    //对动态的各个变量涉及的分布个数，动态离散状态列表，动态条件分布，进行改变
	tdistrnum[k-1]=st3.size();
	tdistrnum[l-1]=st3.size();
	tstate[k-1].clear();
	tstate[l-1].clear();
	for(i=0;i<st3.size();i++){for(j=0;j<st3[i].size();j++){tstate[k-1].push_back(st3[i][j]);tstate[l-1].push_back(st3[i][j]);}}
	tprob[k-1].clear();
	for(i=0;i<prk.size();i++){for(j=0;j<prk[i].size();j++){tprob[k-1].push_back(prk[i][j]);}}
	tprob[l-1].clear();
	for(i=0;i<prl.size();i++){for(j=0;j<prl[i].size();j++){tprob[l-1].push_back(prl[i][j]);}}

	tverdisFa[k-1]=D;                //完成exchange后的tverdisFa[k-1],tverdisFa[k-1]发生了改变
	tverconFa[k-1].clear();
	tverconFa[k-1].push_back(k);
	for(i=0;i<C.size();i++){tverconFa[k-1].push_back(C[i]);} //完成exchange后的tverconFa[k-1]，tverconFa[k-1]发生了改变

	vector<int> Uk(D.size()+C.size());
	merge(D.begin(),D.end(),C.begin(),C.end(),Uk.begin());
	sort(Uk.begin(),Uk.end());
	tverFa[k-1].clear();
    tverFa[k-1].push_back(k);
	for(i=0;i<Uk.size();i++){tverFa[k-1].push_back(Uk[i]);}//完成exchange后的tverFa[k-1],tverFa[k-1]发生了改变
    it=find(tverChi[k-1].begin(),tverChi[k-1].end(),l);
	tverChi[k-1].erase(it);    //在k的子结点集中去掉l,完成exchange后的tverChi[k-1]，tverChi[k-1]发生了改变
	it=tverChi[l-1].begin();   //在l的子结点集中增加k,完成exchange后的tverChi[l-1]，tverChi[l-1]发生了改变
	it++;
	while(it!=tverChi[l-1].end()){if(*it>k){tverChi[l-1].insert(it,k);break;}it++;}
	if(it==tverChi[l-1].end()) tverChi[l-1].push_back(k);
	tverdisFa[l-1]=D;              //完成exchange后的tverdisFa[l-1]，tverdisFa[l-1]发生了改变

	tverconFa[l-1].clear();
	tverconFa[l-1].push_back(l);
	for(i=0;i<CC.size();i++){tverconFa[l-1].push_back(CC[i]);}//完成exchange后的tverconFa[l-1]，tverconFa[l-1]发生了改变
	vector<int> Ul(D.size()+CC.size());
	merge(D.begin(),D.end(),CC.begin(),CC.end(),Ul.begin());
	sort(Ul.begin(),Ul.end());
	tverFa[l-1].clear();
    tverFa[l-1].push_back(l);
	for(i=0;i<Ul.size();i++){tverFa[l-1].push_back(Ul[i]);}//完成exchange后的tverFa[l-1]，tverFa[l-1]发生了改变
	int v;
	for(i=1;i<tverFa[l-1].size();i++){  //exchange后点l的父结点的子结点集发生了改变
		v=tverFa[l-1][i];
		if(find(tverChi[v-1].begin(),tverChi[v-1].end(),l)==tverChi[v-1].end()){
            it=tverChi[v-1].begin();
			it++;
			while(it!=tverChi[v-1].end()){if(*it>l){tverChi[v-1].insert(it,l);break;}it++;}
			if(it==tverChi[v-1].end()) tverChi[v-1].push_back(l);
		}
		if(find(tverChi[v-1].begin(),tverChi[v-1].end(),k)==tverChi[v-1].end()){
			it=tverChi[v-1].begin();
			it++;
			while(it!=tverChi[v-1].end()){if(*it>k){tverChi[v-1].insert(it,k);break;}it++;}
			if(it==tverChi[v-1].end()) tverChi[v-1].push_back(k);
		}
	}

	Vecconm.clear();
    Conm.cverFa=tverFa[k-1];
	Conm.cverFasta=st3;
	Conm.cverFaf=prk;
	Vecconm.push_back(Conm);
	Conm.cverFa=tverFa[l-1];
	Conm.cverFasta=st3;
	Conm.cverFaf=prl;
	Vecconm.push_back(Conm);

    return Vecconm;
}



vector<Dismes> Propagation::Exchange(int k, Dismes & disk, int l, Dismes & disl){
        int i,j;

        vector<Dismes> Vecdism;  //Vecdism[0]描述经过exchange后k的离散信息，Vecdism[1]描述exchange后l的离散信息
        Dismes dism;
        vector<int>::iterator it;
        vector<vector<int> > sk,sl, stk,stl; //stk和stl分别表示exchange后点k和点l涉及的离散状态
        vector<double> pk,pl, prk,prl;//prk和prl分别表示exchange后点k和点l涉及的对应状态下的概率
        vector<int> Fa,lFa,kFa;     //Fa是exchange后k和l共同的父集，都是离散变量，kFa是exchange后k的家族集其中有l
        vector<int> posk,posl,veA,veB,veC;
        double p;
        int posA;

        for(i=1;i<tverFa[k-1].size();i++){Fa.push_back(tverFa[k-1][i]);}
        for(i=1;i<tverFa[l-1].size();i++){if(tverFa[l-1][i]!=k) Fa.push_back(tverFa[l-1][i]);}
        sort(Fa.begin(),Fa.end());
        it=unique(Fa.begin(),Fa.end());
	    Fa.erase(it,Fa.end()); //去掉Fa中多余的点
	    lFa.push_back(l);
	    for(i=0;i<Fa.size();i++){lFa.push_back(Fa[i]);}//exchange后l的家族集
	    Fa.push_back(l);  //在Fa中添加点l
	    sort(Fa.begin(),Fa.end());
	    kFa.push_back(k);
	    for(i=0;i<Fa.size();i++){kFa.push_back(Fa[i]);} //exchange后k的家族集
	    stk=EDisstate(kFa);//exchange后k的家族集对应的状态
	    stl=EDisstate(lFa);//exchange后l的家族集对应的状态

	    sk=disk.dversta;//exchange前k的家族集对应的状态
	    sl=disl.dversta;//exchange前l的家族集对应的状态
	    pk=disk.dverp;//exchange前k的概率分布
	    pl=disl.dverp;//exchange前l的概率分布
        posA=Pos(k,disl.dver);
	    for(i=0;i<stl.size();i++){
            p=0;
            posl=Commpos(lFa,stl[i],disl.dver);  //posl和posk的size()一致,它们的size()应是变量k的取值个数
            posk=Commpos(lFa,stl[i],disk.dver);
            for(j=0;j<posl.size();j++){p=p+pl[posl[j]-1]*pk[posk[j]-1];}
            prl.push_back(p);
        }


        for(i=0;i<stk.size();i++){
            veA=Commpos(kFa,stk[i],disl.dver);
            veB=Commpos(kFa,stk[i],disk.dver);
            veC=Commpos(kFa,stk[i],lFa);


            prk.push_back(pl[veA[0]-1]*pk[veB[0]-1]/prl[veC[0]-1]);
        }
        dism.dver=kFa;
        dism.dversta=stk;
        dism.dverp=prk;
        Vecdism.push_back(dism);
        dism.dver=lFa;
        dism.dversta=stl;
        dism.dverp=prl;
        Vecdism.push_back(dism);


        tdistrnum[k-1]=stk.size();
        tdistrnum[l-1]=stl.size();
	    tstate[k-1].clear();
	    tstate[l-1].clear();
	    for(i=0;i<stk.size();i++){for(j=0;j<stk[i].size();j++){tstate[k-1].push_back(stk[i][j]);}}
	    for(i=0;i<stl.size();i++){for(j=0;j<stl[i].size();j++){tstate[l-1].push_back(stl[i][j]);}}
        tprob[k-1].clear();
	    for(i=0;i<prk.size();i++){tprob[k-1].push_back(prk[i]);}
	    tprob[l-1].clear();
	    for(i=0;i<prl.size();i++){tprob[l-1].push_back(prl[i]);}


	    tverFa[l-1].clear();
	    for(i=0;i<lFa.size();i++){tverFa[l-1].push_back(lFa[i]);}
	    tverdisFa[l-1].clear();
	    for(i=0;i<tverFa[l-1].size();i++){tverdisFa[l-1].push_back(tverFa[l-1][i]);}
	    tverFa[k-1].clear();
	    for(i=0;i<kFa.size();i++){tverFa[k-1].push_back(kFa[i]);}
	    tverdisFa[k-1].clear();
	    for(i=0;i<tverFa[k-1].size();i++){tverdisFa[k-1].push_back(tverFa[k-1][i]);}

	    it=find(tverChi[k-1].begin(),tverChi[k-1].end(),l);
	    tverChi[k-1].erase(it);    //在k的子结点集中去掉l,完成exchange后的tverChi[k-1]，tverChi[k-1]发生了改变
	    it=tverChi[l-1].begin();   //在l的子结点集中增加k,完成exchange后的tverChi[l-1]，tverChi[l-1]发生了改变
	    it++;
	    while(it!=tverChi[l-1].end()){if(*it>k){tverChi[l-1].insert(it,k);break;}it++;}
	    if(it==tverChi[l-1].end()) tverChi[l-1].push_back(k);

	    int v;
	    for(i=1;i<tverFa[l-1].size();i++){  //exchange后点l的父结点的子结点集发生了改变
		  v=tverFa[l-1][i];
		  if(find(tverChi[v-1].begin(),tverChi[v-1].end(),l)==tverChi[v-1].end()){
            it=tverChi[v-1].begin();
			it++;
			while(it!=tverChi[v-1].end()){if(*it>l){tverChi[v-1].insert(it,l);break;}it++;}
			if(it==tverChi[v-1].end()) tverChi[v-1].push_back(l);
		  }
		  if(find(tverChi[v-1].begin(),tverChi[v-1].end(),k)==tverChi[v-1].end()){
			it=tverChi[v-1].begin();
			it++;
			while(it!=tverChi[v-1].end()){if(*it>k){tverChi[v-1].insert(it,k);break;}it++;}
			if(it==tverChi[v-1].end()) tverChi[v-1].push_back(k);
		  }
	    }
        return Vecdism;
}





void Propagation::Messagedis(){
	//cout<<"开始信息分发： "<<endl;
	int i;
	tverFa=verFa;
	tverChi=verChi;
	tverdisFa=verdisFa;
	tverconFa=verconFa;

	tstate=estate;
    tprob=eprob;
	tdistrnum=edistrnum;
	for(i=2;i<cliquenum+1;i++){
//		cout<<"团 "<<i<<" 从它的父团接受信息"<<endl;
	    Dmess(i);
	}
}

void Propagation::Outputstru(int Ver, double cpro, int evidnum, int net_id, std::mt19937& gen) {
	std::ostringstream structureName;
	structureName << "results/structure/structure_" << "V" << Ver << "C" << cpro*100 << "_" << net_id << ".txt";
	std::ofstream outfile(structureName.str());
	if (!outfile) {
		std::cerr << "Unable to save data" << std::endl;
		return;
	}
	else {
		outfile << verFa.size() << endl;
		for (size_t i = 0; i < verFa.size(); i++) {
			if (dcver[verFa[i][0] - 1] == -1) {
				outfile << -1 << " " << 2 << " ";
				for (size_t j = 0; j < verFa[i].size(); j++) { outfile << verFa[i][j] << " "; }outfile << 0 << endl;
			}
			else {
				outfile << -2 << " " << 1 << " ";
				for (size_t j = 0; j < verFa[i].size(); j++) { outfile << verFa[i][j] << " "; }outfile << 0 << endl;
			}
		}

		for (size_t i = 0; i < verFa.size(); i++) {
			outfile << i + 1 << endl; 
			int l = 0; 
			int m = 0;
			for (size_t j = 0; j < distrnum[i]; j++) {
				for (size_t k = 0; k < state[i].size() / distrnum[i]; k++) {
					outfile << state[i][l + k] << " ";
				}
				for (size_t k = 0; k < probcon[i].size() / distrnum[i]; k++) {
					outfile << probcon[i][m + k] << " ";
				}outfile << endl;
				l = l + state[i].size() / distrnum[i];
				m = m + probcon[i].size() / distrnum[i];
			}
		}
	}

	std::ostringstream evidenceName;
	evidenceName << "results/evidence/evidence_" << "V" << Ver << "C" << cpro * 100 << "E" << evidnum << ".csv";
	std::ofstream evidOut(evidenceName.str());
	if (!evidOut) {
		std::cerr << "Unable to save data" << std::endl;
		return;
	}else{
	    // Combine disevid and conevid into evidence_id
		vector<int> evidence_id(disevid);
		evidence_id.insert(evidence_id.end(), conevid.begin(), conevid.end());

		// Combine disevidsta and conevidval into evidence_val
		vector<double> evidence_val;
		for (int val : disevidsta) {
			evidence_val.push_back(static_cast<double>(val));
		}
		evidence_val.insert(evidence_val.end(), conevidval.begin(), conevidval.end());

		// Write combined evidence to the file
		for (size_t i = 0; i < evidence_id.size(); ++i) {
			evidOut << evidence_id[i] << (i < evidence_id.size() - 1 ? "," : "\n");
		}
		for (size_t i = 0; i < evidence_val.size(); ++i) {
			evidOut << evidence_val[i] << (i < evidence_val.size() - 1 ? "," : "\n");
		}
   }
}


void Propagation::OutputNet(int Ver, double cpro, int evidnum, int net_id, std::mt19937& gen) {
	std::ostringstream netName;
	netName << "results/structure/network_" << "V" << Ver << "C" << cpro * 100 << "_" << net_id << ".tgf";
	std::ofstream outnet(netName.str());
	if (!outnet) {
		std::cerr << "Unable to save net file" << std::endl;
		return;
	}
	else {
		for (size_t i = 0; i < verFa.size(); i++) {
			outnet << verFa[i][0] << "\t" << verFa[i][0] << endl;
		}
		outnet << "#" << endl;
		for (size_t i = 0; i < verFa.size(); i++) {
			for (size_t j = 1; j < verFa[i].size(); j++) {
				outnet << verFa[i][j] << "\t" << verFa[i][0] << endl;
			}
		}
	}

	std::ostringstream netName2;
	netName2 << "results/structure/network_" << "V" << Ver << "C" << cpro * 100 << ".m";
	std::ofstream outnet2(netName2.str());
	if (!outnet2) {
		std::cerr << "Unable to save net file for MATLAB" << std::endl;
		return;
	}
	else {
		outnet2 << "our = repmat(struct('filename', {}, 'self', {}, 'parents', {}, 'children', {}, ..." << endl;
		outnet2 << "'values', {}, 'xy', {}, 'discrete', {}, 'index', {}, 'cindex', {}, 'pindex', {}), ..." << endl;
		outnet2 << "1, "<< verFa.size() <<"); " << endl;
		outnet2 << endl;

		for (size_t i = 0; i < verFa.size(); i++) {
			outnet2 << "our(" << i + 1 << ").self = 'V" << i + 1 << "';" << endl;
			for (size_t j = 1; j < verChi[i].size(); j++) {
				outnet2 << "our(" << verChi[i][0] << ").children{" << j << "} = 'V" << verChi[i][j] << "';" << endl;
			}
			for (size_t k = 1; k < verFa[i].size(); k++) {
				outnet2 << "our(" << verFa[i][0] << ").parents{" << k << "} = 'V" << verFa[i][k] << "';" << endl;
			}
			if (dcver[verFa[i][0] - 1] == -1) {
				outnet2 << "our(" << i + 1 << ").discrete = true;" << endl;
				outnet2 << "our(" << i + 1 << ").values = { 1, 2 };" << endl;
			}
			else {
				outnet2 << "our(" << i + 1 << ").discrete = false;" << endl;
			}
			//该属性的作用是指示所在列，与data.csv中的列名一致时可自动读取到，不用设定。
			//outnet2 << "our(" << i + 1 << ").colind = " << i + 1 << ";" << endl;
			outnet2 << endl;
		}
	}
}


void Propagation::OutputData(int Ver, double cpro, int evidnum, int n_sample, int net_id, std::mt19937& gen) {
	std::ostringstream fileName;
	fileName << "results/data/data_" << "V" << Ver << "C" << cpro * 100 << "E" << evidnum << "N" << n_sample << "_" << net_id << ".csv";
	std::ofstream outfile(fileName.str());
	if (!outfile) {
		std::cerr << "Unable to save data" << std::endl;
		return;
	}
	
	// 输出表头
	for (size_t i = 0; i < verFa.size(); ++i) {
		outfile << "V" << (i + 1) << (i < verFa.size() - 1 ? "," : "\n");
	}
	
	// 生成数据
	for (size_t n = 0; n < n_sample; ++n) {
		//cout << "第" << n + 1 << "个样本" << endl;
		std::vector<double> nodesValue;
		for (size_t i = 0; i < verFa.size(); ++i) {
			const auto& node_with_parents = verFa[i];
			vector<int> parentsIndex(node_with_parents.begin() + 1, node_with_parents.end());
			vector<int> contParentsIndex;
			vector<int> parentsState;

			for (size_t idx : parentsIndex) {
				// 根节点没有父节点会跳过
				if (dcver[idx - 1] == -1) { // 只考虑离散父节点
					parentsState.push_back(nodesValue[idx - 1]);
				}
				else {
					contParentsIndex.push_back(idx);
				}
			}

			vector<vector<int>> statesTable;
			vector<vector<double>> parasTable;

			int l = 0;
			int m = 0;
			for (size_t j = 0; j < distrnum[i]; j++) {
				vector<int> possibleStates;
				vector<double> possibleParas;
				for (size_t k = 0; k < state[i].size() / distrnum[i]; k++) {
					possibleStates.push_back(state[i][l + k]);
				}
				for (size_t k = 0; k < probcon[i].size() / distrnum[i]; k++) {
					possibleParas.push_back(probcon[i][m + k]);
				}

				l = l + state[i].size() / distrnum[i];
				m = m + probcon[i].size() / distrnum[i];

				statesTable.push_back(possibleStates);
				parasTable.push_back(possibleParas);
			}
			

			double nodeValue = 100;
			// 离散变量
			if (dcver[node_with_parents[0] - 1] == -1) {
				vector<int> realStates1 = { 1 };
				copy(parentsState.begin(), parentsState.end(), back_inserter(realStates1));

				double realProb1 = 0.0;
				for (size_t t = 0; t < statesTable.size(); t++) {
					if (statesTable[t] == realStates1) {
						realProb1 = parasTable[t][0];
						break;
					}
				}
				double stateProb = RandGen::uniformProb(0, 1, gen);
				if (stateProb <= realProb1) {
					nodeValue = 1;
				}
				else {
					nodeValue = 2;
				}
			}
			// 连续变量
			else {
				vector<double> realParas;
				for (size_t t = 0; t < statesTable.size(); t++) {
					if (statesTable[t] == parentsState) {
						realParas = parasTable[t];
						break;
					}
				}
				
				double mean = realParas.front();
				double sigma2 = realParas.back();

				// 只需要考虑连续父节点
				for (int p = 0; p < contParentsIndex.size(); p++) {
					mean = mean + realParas[p + 1] * nodesValue[parentsIndex[p] - 1]; 
				}
				nodeValue = RandGen::normValue(mean, sqrt(sigma2), gen);
			}
			nodesValue.push_back(nodeValue);
		}

		// 输出数据
		for (size_t i = 0; i < nodesValue.size(); ++i) {
			outfile << nodesValue[i] << (i < nodesValue.size() - 1 ? "," : "\n");
		}
	}
}


Message Propagation::Exchange(int k, Message & mes){
	//cout<<"调用exchange消除点 "<<k<<" "<<endl;
    int i;
	int l,posk,posl;
	Conmes conk,conl;
	Dismes disk,disl;
	vector<Conmes> vecmes=mes.cmes,vecm;
	vector<Dismes> vedmes=mes.dmes,vedm;
	Message mess;

	posk=Posmess(k,mes);
	vector<int>::iterator it;
	vector<int> vec=tverChi[k-1];
    if(dcver[k-1]==-2){
	   for(i=1;i<vec.size();i++){
          conk=vecmes[posk-1];
          l=vec[i];
          posl=Posmess(l,mes);
          conl=vecmes[posl-1];
	      vecm=Exchange(k,conk,l,conl);  //k应尽可能与小的l先exchange
	      vecmes[posk-1]=vecm[0];
          vecmes[posl-1]=vecm[1];
	   }
    }
    else{
         for(i=1;i<vec.size();i++){
             disk=vedmes[posk-1];
             l=vec[i];
             posl=Posmess(l,mes);
             disl=vedmes[posl-1];
             vedm=Exchange(k,disk,l,disl);
             vedmes[posk-1]=vedm[0];
             vedmes[posl-1]=vedm[1];
         }
    }
    for(i=1;i<tverFa[k-1].size();i++){    //点k经过exchange后，成为孤立的点，这样在进一步使用k之外的点时
       l=tverFa[k-1][i];                  //这些点的父结点与子结点中都不会含有k，而k本身记录了完全exchange后，k的父结点
	   it=find(tverChi[l-1].begin(),tverChi[l-1].end(),k);
	   if(it!=tverChi[l-1].end()) tverChi[l-1].erase(it);
	}
	for(i=0;i<vecmes.size();i++){if(vecmes[i].cverFa[0]!=k) mess.cmes.push_back(vecmes[i]);}
	for(i=0;i<vedmes.size();i++){if(vedmes[i].dver[0]!=k) mess.dmes.push_back(vedmes[i]);}

    return mess;
}


void Propagation::Totalmess(){
    //cout<<"调用Totalmess()"<<endl;
    int k;
    int i,j,l;
    vector<int> adjcliq;
    for(k=1;k<cliquenum+1;k++){
        adjcliq=adjcliques[k-1];
        for(i=0;i<adjcliq.size();i++){
            l=adjcliq[i];
            if(l==cliquefa[k-1]){
               if(k!=1){
                 for(j=0;j<smess2[k-1].cmes.size();j++){totalmess[k-1].cmes.push_back(smess2[k-1].cmes[j]);}
                 for(j=0;j<smess2[k-1].dmes.size();j++){totalmess[k-1].dmes.push_back(smess2[k-1].dmes[j]);}
               }
            }
            else{
               for(j=0;j<smess1[l-1].cmes.size();j++){totalmess[k-1].cmes.push_back(smess1[l-1].cmes[j]);}
               for(j=0;j<smess1[l-1].dmes.size();j++){totalmess[k-1].dmes.push_back(smess1[l-1].dmes[j]);}
            }
        }
    }
}


void Propagation::Output(vector<Dispost> & dis, vector<Conpost> & con){
     ofstream outfile("posterior.txt");
     int i,j,k,l;
     if(!outfile) cerr<<"Unable to save result data"<<endl;
     else{
         for(i=0;i<poster.size();i++){
             if(dcver[poster[i]-1]==-1){
                 for(j=0;j<dis.size();j++){
                     if(poster[i]==dis[j].disver){
                        outfile<<"点"<<dis[j].disver<<endl;
                        for(k=0;k<dis[j].state.size();k++){outfile<<dis[j].state[k]<<" ";}outfile<<endl;
                        for(k=0;k<dis[j].dispr.size();k++){outfile<<dis[j].dispr[k]<<" ";}outfile<<endl;
                        break;
                     }
                 }
             }
             else{
                  for(j=0;j<con.size();j++){
                      if(poster[i]==con[j].conver){
                        outfile<<"点"<<con[j].conver<<endl;
                        for(k=0;k<con[j].stap.size();k++){outfile<<con[j].stap[k]<<" ";}outfile<<endl;
                        for(k=0;k<con[j].converf.size();k++){
                            for(l=0;l<con[j].converf[k].size();l++){
                                outfile<<con[j].converf[k][l]<<" ";
                            }outfile<<endl;
                        }outfile<<endl;
                        break;
                      }
                  }

             }
         }

     }

}




void Propagation::Composterior(){
    //cout<<"调用Composterior"<<endl;
    int i,j,k,m;
    vector<int> cliqvex,allvex,rvex,RN,RE; //团中分配的点集,团中所有点的集合,allvex/cliqvex
    vector<int> econ1,econ2;//rvex中消除的连续变量，cliqvex中消除的连续变量
    vector<int> edis1,edis2;//rvex中消除的离散变量，cliqvex中消除的离散变量
    vector<int> never,partall,rem;
    vector<int> elimordering;
    vector<int> relvec;
    vector<vector<int> > disstru;
    Message mess,mess1;
    vector<Conmes> conmes;
	vector<Dismes> dismes,nedism;
	vector<Message> seqmes;
    Dispost dp;  //离散后验
	Conpost cp;  //连续后验
	vector<Dispost> Dpost;//离散点的后验信息
	vector<Conpost> Cpost;//连续点的后验信息
	Totalmess();
    for(i=0;i<cliquefun.size();i++){
       //cout<<"团 "<<i+1<<" 上的性质"<<endl;
       cliqvex=cliquefun[i];
	   if(cliqvex.size()>=1){
             mess=totalmess[i];
             allvex=Allv(mess);
             rvex=Setminus(allvex,cliqvex);
             RN.clear();RE.clear();
             for(j=0;j<rvex.size();j++){if(elabel[rvex[j]-1]==0) RN.push_back(rvex[j]);else RE.push_back(rvex[j]);}
	         Localstru(mess);  //对结构进行重置
             mess=Elimbarren(cliqvex,RN,mess);
             Localstru(mess);
             relvec=Dreach(cliqvex,RE);//点cliqvex能够dreach的变量集
             mess=Exemess(relvec,mess);//得到跟相关变量相关的信息
             Localstru(mess);       //对结构进行重置
             allvex=Allv(mess);
             rvex=Setminus(allvex,cliqvex);
             econ1.clear();econ2.clear();edis1.clear();edis2.clear();
             for(j=0;j<rvex.size();j++){
                 if(dcver[rvex[j]-1]==-2) econ1.push_back(rvex[j]);
                 if((dcver[rvex[j]-1]==-1)&&(elabel[rvex[j]-1]==0)) edis1.push_back(rvex[j]);
			 }
             for(j=0;j<cliqvex.size();j++){
                 if(dcver[cliqvex[j]-1]==-2) econ2.push_back(cliqvex[j]);
                 if((dcver[cliqvex[j]-1]==-1)&&(elabel[cliqvex[j]-1]==0)) edis2.push_back(cliqvex[j]);
			 }
             sort(econ1.begin(),econ1.end());sort(econ2.begin(),econ2.end());
             sort(edis1.begin(),edis1.end());sort(edis2.begin(),edis2.end());
             seqmes.clear();never.clear();
             for(j=econ1.size()-1;j>=0;j--){
                 mess1=mess;
                 if(elabel[econ1[j]-1]==0){mess=Exchange(econ1[j],mess1);}
                 else never.push_back(econ1[j]);
			 }
             seqmes.push_back(mess);
             for(j=econ2.size()-1;j>=0;j--){
                  mess1=mess;
                 if(elabel[econ2[j]-1]==0){mess=Exchange(econ2[j],mess1);seqmes.push_back(mess);}
                 else never.push_back(econ2[j]);
			 }
             mess.cmes.clear();
	         nedism=Dgenbycv(never);
             for(j=0;j<nedism.size();j++){mess.dmes.push_back(nedism[j]);}
             seqmes.push_back(mess);
             if(!edis2.empty()){
                 elimordering=Elimordering(edis2, mess.dmes);
                 for(j=elimordering.size()-1;j>=0;j--){
                     mess1=mess;
                     mess=Elimdisone(elimordering[j],mess1);
                     if(j<edis2.size()+1) seqmes.push_back(mess);
				 }
			 }
             for(j=0;j<cliqvex.size();j++){
                  if(elabel[cliqvex[j]-1]==0){
                     m=0;
                     for(k=0;k<seqmes.size();k++){
                         partall=Allv(seqmes[k]);
                         if(find(partall.begin(),partall.end(),cliqvex[j])!=partall.end()) m++;
                         else break;
					 }
                     if(dcver[cliqvex[j]-1]==-2){cp=Comconp(cliqvex[j],seqmes[m-1]);Cpost.push_back(cp);}
                     else{dp=Comdisp(cliqvex[j],seqmes[m-1]);Dpost.push_back(dp);}
				  }
			 }
	   }
	}
    //Output(Dpost,Cpost);  输出所有的后验分布
}


vector<Dismes> Propagation::Dgenbycv(vector<int> & vec){//返回值是：由条件密度仅是离散变量函数的连续首点集vec产生的离散信息
	//cout<<"调用函数Dgenbycv()"<<endl;
	int i,j,q;
	Dismes cdmes;
	vector<Dismes> Xcdmes;
	double pp;
	vector<double> mp;
	vector<int> pos;
	vector<vector<int> > MB;
    vector<vector<double> > MP;
    for(i=0;i<vec.size();i++){
		cdmes.dver=tverdisFa[vec[i]-1];
		cdmes.dversta=EDisstate(cdmes.dver);
		MB=Tstate(vec[i]);
		MP=Tprob(vec[i]);
		pos.clear();
		for(j=0;j<tverconFa[vec[i]-1].size();j++){pos.push_back(Pos(tverconFa[vec[i]-1][j],conevid));}
        cdmes.dverp.clear();
        for(j=0;j<cdmes.dversta.size();j++){
           for(q=0;q<MB.size();q++){if(equal(MB[q].begin(),MB[q].end(),cdmes.dversta[j].begin())){mp=MP[q];break;}}
           pp=conevidval[pos[0]-1]-mp[0];
		   for(q=1;q<tverconFa[vec[i]-1].size();q++){pp=pp-conevidval[pos[q]-1]*mp[q];}
		   cdmes.dverp.push_back((exp((-1/(2*mp[mp.size()-1]))*pp*pp))/(sqrt(2*3.14*mp[mp.size()-1])));
		}
        Xcdmes.push_back(cdmes);
	}
	//cout<<"调用函数Dgenbycv()完毕"<<endl;
	return Xcdmes;
}


Conpost Propagation::Comconp(int k, Message & mes){
	//cout<<"调用Comconp("<<k<<") "<<endl;
	int i,j;
	Conpost cp;
    vector<int> disfa,confa,allcon,alldis,sta,pos,never;
	double sum=0;
	vector<double> mp,mpp,fp;
	vector<vector<int> > disfasta,alldissta,trucall;
	vector<vector<double> > conprob;
	vector<int> A,R,RN,RE;//A是涉及的全体点，R是A去掉k后的点集，RN是R中非证据变量集,RE是R中证据变量
	A=Allv(mes);R=A;
    vector<int>::iterator it=find(R.begin(),R.end(),k);
    if(it!=R.end()) R.erase(it);
	for(i=0;i<R.size();i++){
        if(elabel[R[i]-1]==1) RE.push_back(R[i]);
        else RN.push_back(R[i]);
    }

    Localstru(mes);
    vector<int> veck,relkvec;
	veck.push_back(k);
    Message remess=Elimbarren(veck,RN,mes),mess,mess1;//remess表示消去barren变量后余下的信息
    Localstru(remess);
	relkvec=Dreach(veck,RE);//点k能够dreach的变量集
    mess=Exemess(relkvec,remess);
    Localstru(mess);

    for(i=0;i<mess.cmes.size();i++){allcon.push_back(mess.cmes[i].cverFa[0]);}
	sort(allcon.begin(),allcon.end());
	for(i=allcon.size()-1;i>=0;i--){
        mess1=mess;
		if(allcon[i]!=k){
			if(elabel[allcon[i]-1]==0) {
                mess=Exchange(allcon[i],mess1);
			}
		    else never.push_back(allcon[i]);
		}
	}
	Exchange(k,mess);

	vector<Dismes> nedism,dism=mess.dmes;
	nedism=Dgenbycv(never);//经过exchange后团中的连续的证据变量的条件密度会成为离散函数



	disfa=tverdisFa[k-1];



    disfasta=EDisstate(disfa);



	conprob=Tprob(k);

	for(i=0;i<nedism.size();i++){dism.push_back(nedism[i]);}
	mess1.dmes.clear();mess1.cmes.clear();
	for(i=0;i<dism.size();i++){mess1.dmes.push_back(dism[i]);}
    vector<int> elimordering=Elimordering(disfa,dism);
    int size=disfa.size();
    for(j=elimordering.size()-1;j>=size;j--){mess1=Elimdisone(elimordering[j],mess1);}//把disfa外的离散变量都加和掉



    dism=mess1.dmes;
	for(i=0;i<disfasta.size();i++){mp.push_back(Multip(disfa,disfasta[i],dism));}
    mpp=mp;
 	for(i=0;i<mpp.size();i++){sum=sum+mpp[i];}
	if(mpp.size()==1) fp.push_back(1);
	else{for(i=0;i<mpp.size();i++){fp.push_back(mpp[i]/sum);}}
    cp.conver=k;cp.stap=fp;



	vector<int> conpos;
	for(i=1;i<tverconFa[k-1].size();i++){conpos.push_back(Pos(tverconFa[k-1][i],conevid));}
	

	
	double mean;
	vector<double> vec;
	vector<vector<double> > converf;
	for(i=0;i<conprob.size();i++){
        vec.clear();
        mean=conprob[i][0];
        for(j=1;j<conprob[i].size()-1;j++){
            mean=mean+conprob[i][j]*conevidval[conpos[j-1]-1]; //??
        }

        vec.push_back(mean);
        vec.push_back(conprob[i][conprob[i].size()-1]);

        converf.push_back(vec);
    }
    cp.converf=converf;
    return cp;
}

vector<int> Propagation::Elimordering(vector<int> & comdisv, vector<Dismes> & dmes){
     vector<vector<int> > disstru=Disstru(dmes);
     Message mess;
     int i,j,k;
     vector<int> elimv;
     for(i=0;i<dmes.size();i++){mess.dmes.push_back(dmes[i]);}
     vector<int> allv=Allv(mess);
     for(i=0;i<allv.size();i++){if(elabel[allv[i]-1]==0) elimv.push_back(allv[i]);}
     vector<int> elimordering;         //顶点集对应的消元顺序
     weight.clear();tweight.clear();label.clear();visit.clear();//weight,tweight,label,visit中顶点顺序与elimv是一致的
     for(i=0;i<elimv.size();i++){weight.push_back(0);}
     for(i=0;i<elimv.size();i++){tweight.push_back(0);}
     for(i=0;i<elimv.size();i++){label.push_back(0);}
     for(i=0;i<elimv.size();i++){visit.push_back(0);}
     vector<int>::iterator iter;
     for(i=0;i<disstru.size();i++){
         if(find(comdisv.begin(),comdisv.end(),disstru[i][0])!=comdisv.end()){
             for(j=0;j<comdisv.size();j++){
                 if(comdisv[j]!=disstru[i][0]){
                    iter=disstru[i].begin();
                    iter++;
                    while(iter!=disstru[i].end()){
                       if(*iter>comdisv[j]){disstru[i].insert(iter,comdisv[j]);break;}
                       iter++;
                    }
                    if(iter==disstru[i].end()) disstru[i].push_back(comdisv[j]);
                 }
             }
         }
     }
     stru=disstru;
     int z,u;
     if(comdisv.empty()){
         for(i=0;i<elimv.size();i++){
             if(i==0){z=elimv[0];label[0]=1;}
             else z=Maxpoint(comdisv,elimv);
             elimordering.push_back(z);
             for(j=0;j<elimv.size();j++){tweight[j]=weight[j];}
             for(j=0;j<elimv.size();j++){
                if(label[j]==0){
                    for(k=0;k<elimv.size();k++){visit[k]=0;}
                    u=Changew(z,elimv[j],elimv);
                    if(u==1) weight[j]++;
                }
             }
         }
     }
     else{
          for(i=0;i<elimv.size();i++){
              if(i==0){z=comdisv[0];label[Pos(z,elimv)-1]=1;}
              else z=Maxpoint(comdisv,elimv);
              elimordering.push_back(z);
              for(j=0;j<elimv.size();j++){tweight[j]=weight[j];}
              for(j=0;j<elimv.size();j++){
                  if(label[j]==0){
                    for(k=0;k<elimv.size();k++){visit[k]=0;}
                    u=Changew(z,elimv[j],elimv);
                    if(u==1) weight[j]++;
                  }
              }

          }
     }
     return elimordering;
}


int Propagation::Maxpoint(vector<int> & comdisv, vector<int> & elimv){
	int i;
	int j=-1,pos,z;
	vector<int> dis;
	if(comdisv.empty()){
	   	for(i=0;i<elimv.size();i++){
		    if((label[i]==0)&&(weight[i]>j)){
                 j=weight[i];
                 pos=i+1;                     //记录权重更大的顶点
            }
	    }
	    z=elimv[pos-1];
	    label[pos-1]=1;
	}
	else{
	   for(i=0;i<elimv.size();i++){
           if((label[i]==0)&&(weight[i]>j)) j=weight[i];
	   }
	   for(i=0;i<elimv.size();i++){
	       if((label[i]==0)&&(weight[i]==j)) dis.push_back(elimv[i]);
	   }
	   for(i=0;i<dis.size();i++){
	       if(find(comdisv.begin(),comdisv.end(),dis[i])!=comdisv.end()){
	          z=dis[i];
	          label[Pos(z,elimv)-1]=1;
	          return z;
           }
	   }
	   z=dis[0];
	   label[Pos(z,elimv)-1]=1;
	}
	return z;
}



int Propagation::Changew(int k, int l, vector<int> & elimv){
    int i=1;
    int m;
    visit[Pos(k,elimv)-1]=1;

    vector<int> vec=stru[Pos(k,elimv)-1];
    while(i!=vec.size()){
        m=vec[i];
        if(m==l) return 1;
        else{
             if((label[Pos(m,elimv)-1]==0)&&(visit[Pos(m,elimv)-1]==0)){
                if(tweight[Pos(m,elimv)-1]<tweight[Pos(l,elimv)-1]){
                   if(Changew(m,l,elimv)==1) return 1;
                }
             }
        }
        i++;
    }
    return 0;
}



void OperonG::Snumbering(){
    int k,l,m,n;         //用来循环
	int z;              //记录最大权的顶点
	int u;
	vector<int>::iterator iter; //用来做泛型迭代
	for(k=0;k<V+1;k++){
        //输出检查
		//cout<<"循环位置 "<<k<<endl;
		if(k==0) {
			z=V+1;
			slabel[V]=1;
			//输出检查
			//cout<<"第一个被编号顶点 "<<V+1<<endl;
		}
		else z=Compare();
		snumbering[z-1]=V+1-k;
		for(l=0;l<V+1;l++){
		    stweight[l]=sweight[l];
		}
		for(m=0;m<V+1;m++){
			if(slabel[m]==0){
				for(n=0;n<V+1;n++){
				      svisit[n]=0;

				}
				u=Updateweight(z,m+1);
				//输出检查，考察星图编号时权重的更新情况
				//cout<<"Updateweight( "<<z<<" , "<<m+1<<" ) "<<u<<endl;
				if(u==1){
                         sweight[m]++;
						 iter=find(striadj[z-1].begin(),striadj[z-1].end(),m+1);
					     if(iter==striadj[z-1].end()){
							 striadj[z-1].push_back(m+1);
							 striadj[m].push_back(z);
						 }

				}
			}

		}

	}
	for(k=0;k<V;k++){
	    mtriadj[k]=striadj[k];
        iter=find(mtriadj[k].begin(),mtriadj[k].end(),V+1);
		if(iter!=mtriadj[k].end())
		    mtriadj[k].erase(iter);
	}
}

int OperonG::Updateweight(int k,int l){    //星图中顶点k被编号的时候，未编号顶点l的权重是否增加
    int i=1;
	int m;
	svisit[k-1]=1;
	vector<int> vec=striadj[k-1];
	while(i!=vec.size()){
		  m=vec[i];
		  if(m==l) return 1;
		  else{
			  if((slabel[m-1]==0)&&(svisit[m-1]==0)){
				  if(stweight[m-1]<stweight[l-1]){
                        if(Updateweight(m,l)==1)
                        return 1;
				  }
			  }

		  }
	     i++;
	}
    return 0;
}


int OperonG::Compare(){
	int i;
	int j=-1,z;
	for(i=0;i<V+1;i++){
		if((slabel[i]==0)&&(sweight[i]>j)){
                 j=sweight[i];
                 z=i+1;                     //记录权重更大的顶点
		}
	}
	slabel[z-1]=1;
    //输出检查
	//cout<<"被编号的顶点 "<<z<<endl;
	return z;
}


vector<vector<int> > Propagation::Disstru(vector<Dismes> & dmes){
       Message mess;
       vector<int> vec1,vec2;
       vector<int>::iterator it;
       vector<vector<int> > disstru;
       int i,j,k;
       for(i=0;i<dmes.size();i++){mess.dmes.push_back(dmes[i]);}
       vector<int> allv=Allv(mess);
       for(i=0;i<allv.size();i++){
          if(elabel[allv[i]-1]==0){
            vec1.clear();
            vec2.clear();
            for(j=0;j<dmes.size();j++){
              if(find(dmes[j].dver.begin(),dmes[j].dver.end(),allv[i])!=dmes[j].dver.end()){
                  for(k=0;k<dmes[j].dver.size();k++){
                      if(elabel[dmes[j].dver[k]-1]==0) vec1.push_back(dmes[j].dver[k]);
                  }
              }
            }
            sort(vec1.begin(),vec1.end());
            it=unique(vec1.begin(),vec1.end());
            vec1.erase(it,vec1.end());
            it=find(vec1.begin(),vec1.end(),allv[i]);
            if(it!=vec1.end()) vec1.erase(it);
            vec2.push_back(allv[i]);
            for(j=0;j<vec1.size();j++){vec2.push_back(vec1[j]);}
            disstru.push_back(vec2);
          }
       }
       return disstru;
}


Dispost Propagation::Comdisp(int k, Message & mes){//计算离散变量的后验分布时，mes中只有离散信息mes.dmes
	//cout<<"调用Comdisp()"<<endl;
	int i,j,m,pos;
	vector<int> allver,ksta,allcon,never;   //allver是mes中所有的离散点,ksta是点k的状态
	vector<vector<int> > allversta;
	Dispost dp;
	double sum=0;
	vector<double> mp,mpp,fp;
    vector<Dismes> dism;
    vector<int> dispk;
    Message mess1=mes;
    dispk.push_back(k);
    vector<int> elimordering=Elimordering(dispk,mes.dmes);
    for(j=elimordering.size()-1;j>=1;j--){mess1=Elimdisone(elimordering[j],mess1);}//把k外的离散变量都加和掉
    dism=mess1.dmes;
//    cout<<"it is too long"<<endl;
    allver=Allv(mess1); //k是allver中唯一的非证据离散变量，所以allversta.size()就是变量k的取值个数
	allversta=EDisstate(allver);
	for(i=0;i<allversta.size();i++){mp.push_back(Multip(allver,allversta[i],dism));} //Multip相当于将Inpos和Multiplep合并在一起
	for(i=0;i<valnum[k-1];i++){ksta.push_back(i+1);mpp.push_back(0);}
	pos=Pos(k,allver);
	for(i=0;i<allversta.size();i++){m=allversta[i][pos-1];mpp[m-1]=mp[i];}
    for(i=0;i<mpp.size();i++){sum=sum+mpp[i];}
    for(i=0;i<mpp.size();i++){fp.push_back(mpp[i]/sum);}
	dp.disver=k;
	dp.state=ksta;
	dp.dispr=fp;
	//cout<<"点 "<<k<<" 的后验概率"<<endl;
	//for(i=0;i<fp.size();i++){cout<<fp[i]<<" ";}cout<<endl;
    return dp;
}




double Propagation::Multiplep(vector<int> inpos, vector<Dismes> & dmes){
    int i;
    double p=1;
    for(i=0;i<inpos.size();i++){p=p*dmes[i].dverp[inpos[i]];}
    return p;

}


Message Propagation::Elimdisone(int k, Message & mess){
      //cout<<"Elimdisone() "<<endl;
      int i,j,m,s;
      double pp;
      vector<int> allver,rever,med,repos;
	  vector<double> mp,remp;
      vector<vector<int> > allversta,reversta;
      vector<int>::iterator it;
      Dismes gedism;      //由于对k求和，新产生的离散信息
      vector<Dismes> redism,hadism; //redism记录余下的离散信息，hadism记录为了消掉k需要处理的离散信息
      for(i=0;i<mess.dmes.size();i++){
          it=find(mess.dmes[i].dver.begin(),mess.dmes[i].dver.end(),k);
          if(it!=mess.dmes[i].dver.end()) hadism.push_back(mess.dmes[i]);
          else redism.push_back(mess.dmes[i]);
      }

      for(i=0;i<hadism.size();i++){for(j=0;j<hadism[i].dver.size();j++){allver.push_back(hadism[i].dver[j]);}}
      sort(allver.begin(),allver.end());
	  it=unique(allver.begin(),allver.end());
	  allver.erase(it,allver.end());
	  it=find(allver.begin(),allver.end(),k);
	  allver.erase(it); //allver代表的是消掉k后的离散函数涉及的变量
      allversta=EDisstate(allver);
      vector<int> sta,posik,powerk,inpos;
      vector<vector<int> > vesta;
      int power;
      for(i=0;i<hadism.size();i++){posik.push_back(Pos(k,hadism[i].dver));}
      for(i=0;i<hadism.size();i++){
          power=1;
          for(j=0;j<hadism[i].dver.size();j++){
              if(hadism[i].dver[j]==k) break;
              if(elabel[hadism[i].dver[j]-1]==0){
                  power=power*valnum[hadism[i].dver[j]-1];
              }
          }
          powerk.push_back(power);
      }
      //cout<<"it is so long"<<endl;
      //cout<<allversta.size()<<endl;
      for(i=0;i<allversta.size();i++){
          pp=0;
          for(s=0;s<valnum[k-1];s++){
             vesta.clear();
             if(s==0){
                for(j=0;j<hadism.size();j++){
                    sta.clear();
                    for(m=0;m<hadism[j].dver.size();m++){
                        if(k==hadism[j].dver[m]) sta.push_back(1);
                        else sta.push_back(allversta[i][Pos(hadism[j].dver[m],allver)-1]);
                    }
                    vesta.push_back(sta);
                }
                inpos=Inpos(vesta,hadism);
                pp=pp+Multiplep(inpos,hadism);
             }
             else{
               for(j=0;j<inpos.size();j++){inpos[j]=inpos[j]+powerk[j];}
               pp=pp+Multiplep(inpos,hadism);
             }
          }
          remp.push_back(pp);
      }
      gedism.dver=allver;
      gedism.dversta=allversta;
      gedism.dverp=remp;
      redism.push_back(gedism);
      mess.dmes=redism;
      return mess;
}







vector<int> Propagation::Allv(Message & mes){
    int i,j;
    vector<int> allv;
    for(i=0;i<mes.cmes.size();i++){
        for(j=0;j<mes.cmes[i].cverFa.size();j++){
            allv.push_back(mes.cmes[i].cverFa[0]);
        }
    }
    for(i=0;i<mes.dmes.size();i++){
        for(j=0;j<mes.dmes[i].dver.size();j++){
            allv.push_back(mes.dmes[i].dver[0]);
        }
    }
    sort(allv.begin(),allv.end());
    vector<int>::iterator iter=unique(allv.begin(),allv.end());
	allv.erase(iter,allv.end());
    return allv;
}






double Propagation::Multip(vector<int> & vec,vector<int> & vecs, vector<Dismes> & dmes){
//	cout<<"调用Multip"<<endl;
	int i,j;
	double p=1;
	vector<int> sta;
	vector<vector<int> > allsta;
	for(i=0;i<dmes.size();i++){
		sta.clear();
		for(j=0;j<dmes[i].dver.size();j++){sta.push_back(vecs[Pos(dmes[i].dver[j],vec)-1]);}
		allsta.push_back(sta);
	}
	int pos,power;
	vector<int> inpos;         //利用十进制和多进制的转换关系
	for(i=0;i<allsta.size();i++){
	    pos=0;power=1;
	    for(j=0;j<allsta[i].size();j++){
	        if(elabel[dmes[i].dver[j]-1]==0){
	            pos=pos+(allsta[i][j]-1)*power;
	            power=power*valnum[dmes[i].dver[j]-1];
	        }
	    }
	    inpos.push_back(pos);
	}
	//for(i=0;i<inpos.size();i++){cout<<inpos[i]<<" ";}cout<<endl;
	for(i=0;i<inpos.size();i++){p=p*dmes[i].dverp[inpos[i]];}
//	cout<<"完成Multip"<<endl;
    return p;
}



vector<int> Propagation::Sep(int k,int l){
	int i;
	int ms;
	vector<int> S;
	if(cliquefa[k-1]==l){ms=msignver[k-1];for(i=1;i<mmadj[ms-1].size();i++){S.push_back(mmadj[ms-1][i]);}}
	else{ms=msignver[l-1];for(i=1;i<mmadj[ms-1].size();i++){S.push_back(mmadj[ms-1][i]);}}
    return S;
}

int Propagation::Posmess(int k, Message & mes){
	int i,m=0;
	if(dcver[k-1]==-2){
	    for(i=0;i<mes.cmes.size();i++){
	    m++;
		if(mes.cmes[i].cverFa[0]==k) break;
	   }
    }
	else{
         for(i=0;i<mes.dmes.size();i++){
           m++;
           if(mes.dmes[i].dver[0]==k) break;
        }
    }
    return m;
}



int Propagation::Pos(int k,vector<int> & vec){
    int i,m;
	m=0;
	for(i=0;i<vec.size();i++){m++;if(vec[i]==k) break;}
	return m;
}



Propagation::~Propagation(){
  cout<<"the destructor of Propagation "<<endl;
}



vector<int> Propagation::Remainpos(int k, vector<vector<int> > & vecs){
    int i,j,m;
	vector<int> vecpos,evec,pos,epos,vec;
	vec=tverdisFa[k-1];
	for(i=0;i<vec.size();i++){
		if(elabel[vec[i]-1]==1){
			   evec.push_back(vec[i]);
			   pos.push_back(i+1);
			}
	}
	for(i=0;i<evec.size();i++){epos.push_back(Pos(evec[i],disevid));}
	for(i=0;i<vecs.size();i++){
		m=0;
		for(j=0;j<evec.size();j++){if(disevidsta[epos[j]-1]!=vecs[i][pos[j]-1]){m++;break;}}
		if(m==0){vecpos.push_back(i+1);}
	}
    return vecpos;
}

vector<vector<int> > Propagation::Remainsta(int k){
     int i;
	 vector<int> vecpos;
     vector<vector<int> > sta,rsta;
	 sta=Tstate(k);
	 vecpos=Remainpos(k,sta);
	 for(i=0;i<vecpos.size();i++){rsta.push_back(sta[vecpos[i]-1]);}
	 return rsta;
}

vector<vector<double> > Propagation::Remainpr(int k){
     int i;
	 vector<int> vecpos;
	 vector<vector<int> > sta;
	 vector<vector<double> > pr,rpr;
	 sta=Tstate(k);
	 vecpos=Remainpos(k,sta);
	 pr=Tprob(k);
	 for(i=0;i<vecpos.size();i++){rpr.push_back(pr[vecpos[i]-1]);}
	 return rpr;
}

vector<vector<int> > Propagation::Tstate(int k){
	int i,j,l;
	vector<int> vec,vecpos;
	vector<vector<int> > ve,tsta;
	for(i=0;i<tdistrnum[k-1];i++){
		vec.clear();
		for(j=0;j<tstate[k-1].size()/tdistrnum[k-1];j++){
			l=i*(tstate[k-1].size()/tdistrnum[k-1])+j;
	        vec.push_back(tstate[k-1][l]);
		}

		ve.push_back(vec); //因为这是一个插入函数，即使当tstate[k-1].size()=0时，仍有ve.size()=1
	}
	return ve;
}



vector<vector<double> > Propagation::Tprob(int k){
	int i,j,l;
	vector<double> pvec;
	vector<vector<double> > pve;
	for(i=0;i<tdistrnum[k-1];i++){
	    pvec.clear();
		for(j=0;j<tprob[k-1].size()/tdistrnum[k-1];j++){
		    l=i*(tprob[k-1].size()/tdistrnum[k-1])+j;
			pvec.push_back(tprob[k-1][l]);
		}
		pve.push_back(pvec);
	}
	return pve;
}




int Propagation::Tpos(int k,int l){
        int i;
		for(i=0;i<tverconFa[l-1].size();i++){
		    if(tverconFa[l-1][i]==k)
				 break;
		}
		return i+1;
}








int OperonG::Complete(vector<int> vec){   //检验点集在道义图的极小M三角化图中是否完全
    int l,m;
	int n,s;
	vector<int>::iterator iter1,iter2;
	for(l=0;l<vec.size();l++){
		for(m=l+1;m<vec.size();m++){
	              n=vec[l];
				  s=vec[m];
			      iter1=find(mtriadj[n-1].begin(),mtriadj[n-1].end(),s);
				  iter2=find(mtriadj[s-1].begin(),mtriadj[s-1].end(),n);
				  if((iter1==mtriadj[n-1].end())||(iter2==mtriadj[s-1].end()))
					  return 0;
		}
	}
    return 1;
}

int OperonG::Mverpos(){  //记录下点与团的位置的关系,并返回道义图的极小m三角化图中团的个数
	 int k;
	 int z;
	 int l=0,m=0;      //m最后将表示M素块的总数
	 msignver.clear();
	 mcliquever.clear();
	 for(k=0;k<V;k++){
		 z=morderver[k];     //z是编号为k+1的点
		 //输出检查
		 //cout<<"z "<<z<<" "<<" m "<<m<<" mmadj[z-1] "<<mmadj[z-1].size()<<endl;
         if(l<=mmadj[z-1].size()){
			  m++;
	          mverpos[z-1]=m;             //记录点z所在的团的排号
			  mcliquever.push_front(z);   //记录新团开始的点
			  if(l!=0) msignver.push_front(morderver[k-1]); //记录前一个团的标识点
		 }
		 else mverpos[z-1]=m;
         l=mmadj[z-1].size();
	 }
	 msignver.push_front(morderver[V-1]);
     for(k=0;k<V;k++){
	     mverpos[k]=m+1-mverpos[k];         //点k+1所在的团的排号
	 }
	 //输出检验
/*     for(k=0;k<m;k++){
	   cout<<msignver[k]<<" ";
	 }*/
	 return m;
}

void OperonG::Constructtree(){         //构建道义图的极小m三角化图中团的junction tree
    cliques.clear();
    cliquefa.clear();
    adjcliques.clear();
    cliquenum=Mverpos();
	int k,l;
	int z,x;
    vector<int> vec;
	cliquefa.push_back(0);
	vector<int>::iterator iter;
	for(k=1;k<cliquenum;k++){
		vec.clear();
        z=msignver[k];
        iter=mmadj[z-1].begin();
		for(;iter!=mmadj[z-1].end();iter++){
		    if(iter!=mmadj[z-1].begin())
			   vec.push_back(mnumbering[*iter-1]);    //vec记录当时的单调邻集中点的编号
		}
        iter=min_element(vec.begin(),vec.end());       //取出最小的编号
        x=morderver[*iter-1];                          //最小编号对应的点
		cliquefa.push_back(mverpos[x-1]);                //存入最小编号对应的点所在团的排号
	}
    for(k=0;k<cliquenum;k++){
	    z=mcliquever[k];
	    cliques.push_back(mmadj[z-1]);
	}
    for(k=0;k<cliquenum;k++){
	   vec.clear();
	   vec.push_back(cliquefa[k]);
	   for(l=0;l<cliquenum;l++){
		      if(cliquefa[l]==(k+1)) vec.push_back(l+1);
	   }
	   adjcliques.push_back(vec);
	}
}

void OperonG::Morsmcs(){
	int k,l;//做循环用
	int z;         //记录最大权的顶点
	mmadj.clear();

	for(k=0;k<V;k++){
	   if(k==0) {
		   for(l=0;l<V;l++){
		      if(dcver[l]==-1) break;
		   }
		   if(l!=V){z=l+1; mlabel[z-1]=1;}
		   else{z=1;mlabel[z-1]=1;}
	  //   cout<<"第一个被编号的顶点 "<<z<<endl;
	   }  //如果有离散点，第一个编号点的位置是任一离散点，如果没有，选第一个连续点最先编号
	   else z=Mcompare();
	   mnumbering[z-1]=V-k;
	   morderver.push_front(z);
	   Mupdate(z);
	}
	//输出检查
	int m;
	vector<int> vec;
	for(k=0;k<V;k++){                       //构建单调邻集
		   vec.clear();
	       vector<int>::iterator iter=mtriadj[k].begin();
		   for(;iter!=mtriadj[k].end();iter++){
			    m=*iter;
			    if(mnumbering[m-1]>=mnumbering[k])
		            vec.push_back(m);
		   }
		   mmadj.push_back(vec);
	}

	//输出检查，输出单调邻集
/*	cout<<"单调邻集: "<<endl;
	 for(k=0;k<V;k++){
	    vector<int>::iterator it=mmadj[k].begin();
	   for(;it!=mmadj[k].end();it++){
	        cout<<*it<<" ";
	   }
	   cout<<endl;
	}*/
}

void OperonG::Mupdate(int k){   //道义图中顶点k被编号的时候，增加其他顶点的权重
    int i=1,j;
	vector<int> vec=mtriadj[k-1];

	while(i!=vec.size()){
		 j=vec[i];
		 if(mlabel[j-1]==0) mweight[j-1]++;
	     i++;
	}
}

int OperonG::Mcompare(){
    int i;
	int j=-1,k,z;
	for(i=0;i<V;i++){
		if((mlabel[i]==0)&&(mweight[i]>j)){
                 j=mweight[i];      //j记录的是未编号顶点中第一个最大的权重
                 z=i+1;             //z记录的是第一个权重最大的未编号顶点
		}
	}
	for(i=z-1;i<V;i++){
        k=mweight[i];
		if((k==j)&&(dcver[i]==-1)){
			if(mlabel[i]==0){
			//	   cout<<"被编号的离散顶点 "<<i+1<<endl;
                   mlabel[i]=1;
				   return i+1;
			}
		}
	}
//	cout<<"被编号的连续顶点 "<<z<<endl;
	mlabel[z-1]=1;
	return z;

}


void StarMG::inStarMG() {
	staradj.clear();
	int k;
    for(k=0;k<V;k++){
	  staradj.push_back(moradj[k]);
	}
	vector<int> vec;
	vec.push_back(V+1);
	for(k=0;k<V;k++){
		if(dcver[k]==-1){
		    staradj[k].push_back(V+1);
			vec.push_back(k+1);
		}
	}
	staradj.push_back(vec);
}




void MoralG::inMoralG(){
	int k,l,m;
	int pos1,pos2;
	vector<int>::iterator iter1,iter2;
	moradj.clear();
	for(k=0;k<V;k++){
	   moradj.push_back(verFa[k]);
	}

	vector<int> ver;
	for(k=0;k<V;k++){
	      ver=verFa[k];
		  for(l=0;l<ver.size();l++){
			  for(m=l+1;m<ver.size();m++){     //把家族完全化
			      pos1=ver[l];
				  pos2=ver[m];
				  iter1=find(moradj[pos1-1].begin(),moradj[pos1-1].end(),pos2);
				  iter2=find(moradj[pos2-1].begin(),moradj[pos2-1].end(),pos1);
				  if(iter1==moradj[pos1-1].end()) moradj[pos1-1].push_back(pos2);
				  if(iter2==moradj[pos2-1].end()) moradj[pos2-1].push_back(pos1);
			  }
		  }
	}

}


bool ClearDir(const std::wstring& path) {
	WIN32_FIND_DATA findData;
	HANDLE hFind = FindFirstFile((path + L"\\*").c_str(), &findData);

	if (hFind == INVALID_HANDLE_VALUE) return false;

	do {
		std::wstring fileName = findData.cFileName;
		if (fileName == L"." || fileName == L"..") continue;
		std::wstring fullPath = path + L"\\" + fileName;

		if (findData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {
			ClearDir(fullPath); // Recursively clear directory
			RemoveDirectory(fullPath.c_str()); // Delete the directory itself
		}
		else {
			DeleteFile(fullPath.c_str()); // Delete files directly
		}
	} while (FindNextFile(hFind, &findData) != 0);

	FindClose(hFind);
	return true;
}

bool CreateDir(const std::wstring& path) {
	if (!CreateDirectory(path.c_str(), NULL)) {
		if (GetLastError() == ERROR_ALREADY_EXISTS) {
			return ClearDir(path);
		}
		return false;
	}
	return true;
}

void CreateSubDirs(const std::wstring& base, const std::vector<std::wstring>& dirs) {
	if (!base.empty() && !CreateDir(base)) {
		std::wcout << L"Failed to create or clear base dir: " << base << std::endl;
		return;
	}

	for (const auto& dir : dirs) {
		std::wstring fullPath = base.empty() ? dir : base + L"\\" + dir;
		if (CreateDir(fullPath)) {
			std::wcout << L"Created or cleared: " << fullPath << std::endl;
		}
		else {
			std::wcout << L"Failed: " << fullPath << std::endl;
		}
	}
}


int main() {
	// 设置输出的目录
	std::wstring base = L"results";
	std::vector<std::wstring> dirs = { L"structure", L"data", L"evidence"};
	CreateSubDirs(base, dirs);

	// 模拟设置
	vector<int> samples = { 1, 1000 };
	vector<int> vertexes = { 50, 75, 100 };
	vector<double> continuous = { 0, 0.25, 0.5, 0.75, 1 };
	int n_evidence;
	clock_t time_start, time_mid, time_end;
	ofstream time_log("results/timecost.txt");
	if (!time_log) {
		cerr << "Unable to save time data" << endl;
		return 1;
	}

	Propagation bt;
	for (int n_vertex : vertexes) {
		for (double p_continuous : continuous) {
			// Iterate over each evidence variable count.
			for (size_t e = 0; e < 10; ++e) {
				n_evidence = (e * n_vertex) / 10;
				time_log << n_vertex << " " << p_continuous << " " << n_evidence << endl;
				cout << "顶点数: " << n_vertex << " 连续变量比例: " << p_continuous << " 证据变量数: " << n_evidence << endl;
				// Randomly generate 10 Bayesian networks.
				for (size_t b = 0; b < 100; ++b) {
					std::mt19937 gen(b); // 以b为种子的Mersenne Twister生成器
					bt.inBayesnet(n_vertex, p_continuous, n_evidence, gen);
					bt.Outputstru(n_vertex, p_continuous, n_evidence, b, gen);
					bt.OutputNet(n_vertex, p_continuous, n_evidence, b, gen);
					for (int n_sample : samples) {
						bt.OutputData(n_vertex, p_continuous, n_evidence, n_sample, b, gen);
					}

					time_start = clock();
					bt.inMoralG();
					bt.inStarMG();
					bt.inOperonG();
					time_mid = clock();

					bt.inPropagation();
					bt.Messagecol();
					bt.Messagedis();
					bt.Composterior();
					time_end = clock();

					time_log << static_cast<double>(time_mid - time_start) / CLOCKS_PER_SEC << " "
						<< static_cast<double>(time_end - time_mid) / CLOCKS_PER_SEC << endl;
				}
			}
		}
	}

	system("pause");
	return 0;
}
