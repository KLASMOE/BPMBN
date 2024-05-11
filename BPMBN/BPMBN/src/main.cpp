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
	int expectedEdges = static_cast<int>(numVertices * log(log(numVertices))); // ��������

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
	//���ݶ��������������������Լ�֤�ݱ����������Ա�Ҷ˹���еĸ���Ԫ��дֵ
    void inBayesnet(int Ver, double cpro, int evidnum, std::mt19937& gen);
	void Print();     //��������ͼ
protected:
   vector<vector<int> > verFa;    //V�����㼰�丸�ڵ㼯��ÿһ���Ǳ�Ҷ˹���Ļ����ṹ��������,������С�ŵ���Ŵ洢
   vector<vector<int> > verChi;   //verChi[k]��ʾ����k+1�����ӽ��ļ���
   vector<vector<int> > verdisFa; //V������ļ�������ɢ������
   vector<vector<int> > verconFa; //V������ļ���������������
   vector<vector<int> > state;    //V�����㼰�丸�ڵ㼯�о������ɢ״̬�б�,����������ɢ������ȡ˳���й�
   vector<vector<double> > prob;  //V������ֲ���Ϣ�б�,�����������������Ĵ�ȡ˳���й�
   vector<vector<double> > probcon;

   vector<int> distrnum;  //��������ӵ�еķֲ��������ɼ�������ɢ����ȡ����ֵ����
   vector<int> dcver;     //����������������ɢ�ģ�����������
   vector<int> valnum;    //������ɢ������ȡֵ����������������ȡֵ������Ϊ1

   	vector<int> elabel;           //֤�ݱ�ǩ��elabel[k]=0��ʾ����k+1����֤�ݣ�elabel[k]=1��ʾ����k+1��֤��
    vector<int> disevid;          //��ɢ֤�ݱ�����
	vector<int> disevidsta;       //��ɢ֤�ݱ�����״̬
	vector<int> conevid;          //����֤�ݱ�����
	vector<double> conevidval;    //����֤�ݱ�����ȡֵ


   int V;                //�������
};


void Bayesnet::inBayesnet(int Ver, double cpro, int evidnum, std::mt19937& gen){
       //cout<<"������Ϊ "<<Ver<<" ������������Ϊ "<<cpro<<" ֤�ݱ�������Ϊ "<<evidnum<<" �������Ҷ˹��:"<<endl;
       int i,j,m,n,k,l;
       vector<int> vec,ve;
       vector<double> vec1,vec2,vec3;
       V=Ver;
	   verFa = GenerateRandomDAG(V, gen);
       int conv=static_cast<int> (V*cpro);  //������������
       int disv=V-conv;                     //��ɢ��������,����ɢ��������ǰ��

       verChi.clear();
       verdisFa.clear();
       verconFa.clear();
       state.clear();
       prob.clear();
	   probcon.clear();

       distrnum.clear();
       dcver.clear();
       valnum.clear();
       //�ñ��С�ı�������ɢ��������Ŵ�ı�������������
       for(i=0;i<V;i++){if(i<disv) dcver.push_back(-1);else dcver.push_back(-2);}
       //��������ɢ������ȡֵ������Ϊ2
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
		   if (dcver[i] == -1) { //��ɢ�������ʵĲ�������
			   vec2.clear();
			   for (j = 0; j < distrnum[i]; j++) {
				   vec2.push_back(RandGen::uniformProb(0, 1, gen));
			   }
			   acc = accumulate(vec2.begin(), vec2.end(), 0.0); //�����ۻ��ܺ�

			   // ʹ�� std::transform ��� vec1
			   vec1.reserve(vec2.size()); // Ԥ��������ռ�
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
		   else {  //�����������ʵĲ�������
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
    for(i=0;i<V;i++){tv.push_back(i+1); }//���ö��㼯
	std::shuffle(tv.begin(), tv.end(), gen); //�������֤�ݱ�����

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


class MoralG : public Bayesnet{     //����ͼMoralG,�Ǵ�Bayesnet�̳ж���
public:
       MoralG(){cout<<"the constructor of MoralG()"<<endl;};
	   void inMoralG();
	   ~MoralG(){cout<<"the deconstructor of MoralG()"<<endl;};
protected:
       vector<vector<int> > moradj; //moradj[k]��ʾ����ͼ�ж���k+1����k+1���ڵĶ��㼯
};

class StarMG : public MoralG{       //��ͼ�Ǵӵ���ͼ�̳ж���
public:
       StarMG(){cout<<"the constructor of StarMG()"<<endl;};
       void inStarMG();
	   ~StarMG(){cout<<"the deconstructor of StarMG()"<<endl;};
protected:
	   vector<vector<int> > staradj;//�����ǽڵ�V+1,staradj[V]���ǽڵ㼰�����ڵĶ���
};

class OperonG: public StarMG{       //ͼ�ϵĲ���,������������ͼ�ļ�Сm���ǻ�ͼ���ŵ�junction tree
public:
	 OperonG(){cout<<"the constructor of OperonG()"<<endl;};
	 void inOperonG();
	 void Snumbering();                 //����ͼʹ��mcs-m�㷨����
     int Compare();                     //�Ƚ�Ȩ�ش�С������Ȩ�����ĵ�
	 int Updateweight(int k, int l);    //��ͼ�ж���k����ŵ�ʱ��δ��Ŷ���l��Ȩ���Ƿ�����

	 void Morsmcs();           //�Ե���ͼ�ļ�СM���ǻ�ͼʹ��Smcs�㷨
     int Mcompare();           //�Ƚ�Ȩ�ش�С������Ȩ�����ĵ㣬�����ɢ����������Ȩ��һ����󣬷�����ɢ��
	 void Mupdate(int k);      //����k�����ʱ����������δ��ŵ��ڵ��Ȩ��

     int Complete(vector<int> vec);      //����㼯�ڵ���ͼ�ļ�СM���ǻ�ͼ���Ƿ���ȫ

	 void Constructtree();        //��������ͼ�ļ�Сm���ǻ�ͼ���ŵ�junction tree
	 int Mverpos();               //��¼�µ����ŵ�λ�õĹ�ϵ,�����ص���ͼ�ļ�Сm���ǻ�ͼ���ŵĸ���

	 ~OperonG(){cout<<" the destructor of OperonG: "<<endl;};
protected:
	 vector<int> sweight;           //��ͼ�ж���ĵ�ǰȨ��
	 vector<int> stweight;          //��ͼ�ж������ʱȨ��
	 vector<int> slabel;            //��ͼ�ж�����(1)��(0)�Ѿ������
	 vector<int> snumbering;        //��ͼ�ж���ı�ţ�snumbering[k]��ʾ����k+1�ı��
	 vector<int> svisit;            //��ͼ�ж�����(1)��(0)�����ʹ�
	 vector<vector<int> > striadj;  //striadj[k]��ʾ����ͼ�ļ�С���ǻ�ͼ�е�k+1��k+1�����ڵ㼯

     vector<vector<int> > mtriadj;  //mtriadj[k]��ʾ�ڵ���ͼ�ļ�СM���ǻ�ͼ�е�k+1��k+1�����ڵ㼯
	 vector<vector<int> > mmadj;    //mmadj[k]��ʾ����ͼ�ļ�Сm���ǻ�ͼ�е�k+1��k+1�ĵ����ڼ�
	 vector<int> mnumbering;        //mnumbering[k]��ʾ��k+1�ı��
	 deque<int> morderver;          //��˳���ŵĶ�����, morderver[k]��ʾ���Ϊk+1�ĵ�
	 vector<int> mlabel;            //mlabel[k]��ʾ����k+1��(1)��(0)�Ѿ������
	 vector<int> mweight;           //����ͼ�ж���ĵ�ǰȨ��


	 deque<int> msignver;            //msignver[k]��ʾһ���㣬������ʶ���ź�Ϊk+1���ţ������ĵ����ڼ��Ƿֽ���
	 deque<int> mcliquever;          //mcliquever[k]��ʾһ���㣬����㼰�䵥���ڼ��γ��ź�Ϊk+1����
     vector<vector<int> > cliques;   //�ŵ��б�,cliques[k]��ʾ�ź�Ϊk+1����

	 vector<int> mverpos;            //��ʾʹ���㷨SMCS����D-numbering�¸����㱻��������Ǹ�Ψһ���ŵ��źţ�mverpos[k]��ʾ��k+1�����ŵ��ź�
	 vector<int> cliquefa;           //cliquefa[k]��ʾһ���źţ������������ź�Ϊk+1���ŵĸ��ŵ��ź�
	 vector<vector<int> > adjcliques;//adjcliques[k]��ʾһ���źż����������������ź�Ϊk+1���ڵ������ŵ��ź�
	 int cliquenum;                  //�ŵĸ���
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
	vector<int> dver;            //��ɢ�����ļ��壬��Щ����������֤�ݱ���,����Composterior()��dver������ʾһ����ɢ������
	vector<vector<int> > dversta;//dversta[k]��ʾ��ɢ������dver��k+1����ɢ״̬,��ȡ��dversta��֤�ݱ���״̬�Ժ��״̬�ռ�
	vector<double> dverp;        //dverp[k]��ʾ�����ڵ�k+1����ɢ״̬�µ�ȡֵ,����Composterior()��dverp������ʾһ��ʵ����
};

class Conmes{
public:
	vector<int> cverFa;              //���������ļ���
    vector<vector<int> > cverFasta;  //���������ĸ�������ɢ������״̬�б�,��ȡ����ɢ֤�ݱ���״̬�Ժ��״̬�б�
	vector<vector<double> > cverFaf; //��Ӧ��ɢ״̬�б�������ֲ��о�ֵϵ���������б�
};

class Message{
public:
	vector<Conmes> cmes; //��֮�䴫�ݵ�������������Ϣ,��¼���������ܶȺ�����Ϣ��ÿ�������ܶȺ�������һ����Ϣ
    vector<Dismes> dmes; //��֮�䴫�ݵ���ɢ��������Ϣ����¼�����������ʺ�����Ϣ��ÿ���������ʺ�������һ����Ϣ
};


class Dispost{
public:
	int disver;              //����
	vector<int> state;    //state[k]��ʾ����ĵ�k+1��״̬
	vector<double> dispr;    //pr[k]��ʾ�����ڵ�k+1��״̬�µĺ������
};

class Conpost{
public:
    int conver;                     //��������
	vector<double> stap;            //�����ĺ���ֲ������ɸ�״̬�µĻ�Ϸֲ���stap[k]��ʾ�����ڵ�k+1��״̬�µĸ���
    vector<vector<double> > converf;  //conf[k]��ʾ�����ڵ�k+1��״̬�µĺ����ֵ�ͷ��
};


class Propagation: public OperonG{
public:
    Propagation(){cout<<"the constructor of Propagation()"<<endl;};
	void inPropagation();

    void Messagecol();                             //��Ϣ�ռ�����
	Message Cmess(int k);                          //����ֵ�ǣ��ռ���Ϣʱ����k�����ŵ���Ϣpotential
	Message Proj(Message & mes, int k, int l);     //��һ����k����ϢͶӰ����l�ϣ�����ֵҲ����Ϣpotential,mes��ʾ��k�ϵĵ�ǰ��Ϣpotential
	vector<int> Sep(int k,int l);                  //����ֵ�ǣ���k����l�ķ�����
    Message Projmes(Message & mes, int k, int l); //����ֵ�ǣ����ݵĲ�����Ϣ,�����������ʹ��exchange��ȥ����
//    Message Elimbarren(vector<int> & RN, Message & mes); //��ȥ��֤�ݱ�����RN�е�barren���������ڽṹ�й�����Щ����,ʹ���������ǰ��Ҫ��mes���ṹ����
    Message Elimbarren(vector<int> & TAR, vector<int> & RN, Message & mes);

    //����ֵ�ǣ��ӵ㼯S�����ܹ�dreach�Ķ��㼯(������֤�ݱ���)��vec�Ǹ����ľֲ�֤�ݱ�������
    vector<int> Dreach(vector<int> & S, vector<int> & vec);
    //vec�Ǹ����ľֲ�֤�ݱ���������l��һ����ײ�㣬�Ƿ����l��һ������d��ͨ������ֵ��1��ʾ���ܣ�����ֵ��0���ƺ������ܡ�
    int Acceptl(int l, vector<int> & vec);
    //����l�ĺ�������Ƿ���vec�е��ཻ��1��ʾ�ཻ��0��ʾ���ཻ��
    vector<int> Des(int k);
    //vec�Ǹ����ľֲ�֤�ݱ�����,���ӽ�㵽�����ı�kl�������¼�󣬽�һ�������ж�l���丸���ӽ���γɵı��Ƿ�������¼��
    void Dreachf(int k, int l, vector<int> & vec);
    //vec�Ǹ����ľֲ�֤�ݱ�����,�ɸ���㵽�ӽ��ı�kl�������¼�󣬽�һ�������ж�l���丸�ӽ���γɵı��Ƿ�������¼��
    void Dreachc(int k, int l, vector<int> & vec);
    //����d���������exchange����Ϣ,dsep�Ǳ�d������Ĳ��֣���mes��ȥ��������desp�йص���Ϣ��
    Message Exemess(vector<int> & rels, Message & mes); //����mes��ĳ��potential�漰�ı�������rels�ཻ�ǿգ��������potential


	int Pos(int k, vector<int> & vec);           //����ֵ�ǣ�����k��vec�е�λ��


	vector<int> Remainpos(int k, vector<vector<int> > & vecs);
	//���ر���k�ļ�������ɢ��������״̬vecs�£�������ɢ֤��ȡֵ֮������ɢ֤�ݱ�����һ�µ�״̬��vecs�е�λ��
	vector<vector<int> > Remainsta(int k);
	//�����Ա���kΪ�׵�������ܶȣ��ڸ�����ɢ֤�ݱ���ȡֵ���漰����ɢ״̬
	vector<vector<double> > Remainpr(int k);
	// �����Ա���kΪ�׵�������ܶȣ��ڸ�����ɢ֤�ݱ���ȡֵ���漰����ɢ״̬��Ӧ���ܶȺ���


	vector<vector<int> > Tstate(int k);   //���ص�k�Ķ�̬��ɢ״̬�б����Ǹ����е���ɢ֤�ݱ����Ѿ�������
	vector<vector<double> > Tprob(int k); //���ص�k�Ķ�̬�����ֲ��б����Ǹ����е���ɢ֤�ݱ����Ѿ�������
    int Tpos(int k,int l);  //��������k,l,k��l�ĸ����l��������k��tverconFa[l-1]�е�λ��
	int Posmess(int k, Message & mes); //��k������������������������k����Ϣmes��������Ϣmes.cmes�е�λ��
	//��k����ɢ������������ɢ��������ɢ��Ϣmes.dmes�е�λ��


	vector<Conmes> Exchange(int k, Conmes & conk, int l, Conmes & conl);  //����exchange��ĵ�k�͵�l��������Ϣ
	//�ں�������״̬�У�����������k������tverChi����������С����������l����exchange
    vector<Dismes> Exchange(int k, Dismes & disk, int l, Dismes & disl);  //����exchange��ĵ�k�͵�l����ɢ��Ϣ
    //�ں�������״̬�У�����ɢ����k������tverChi����������С����ɢ����l����exchange
    vector<int> Commonpos(vector<int> & vec, vector<int> & vec_val,
                          vector<int> & ved, vector<vector<int> > & vedsta); //����ɢ������vedΪ����
    //���ر�����vec��ȡֵvec_val��vec��ved�Ľ����Ϲ�ͬ��״̬��ved��״̬�е�λ�ã�����vecsta��vedsta
    //Ϊvec��ved��Ӧ��״̬��
    vector<int> Commpos(vector<int> & vec, vector<int> & vec_val,  vector<int> & ved);

 	Message Exchange(int k, Message & mes);//�ں����϶Ա���kʹ��exchange����,��������������Ϣ,ʹ���������ǰ��Ҫ��mes���ṹ����
    vector<vector<int> > EDisstate(vector<int> & vec);//������ɢ������vec��֤�ݱ���ȡ�������ɢ״̬�б�
    void Localstru(Message & mes);   //�ֲ��ṹ��������ͨ����Ϣmes����дmes�漰�ı����ľֲ��ṹ
    int Domnum(int k);        //���أ�Ϊ��exchange������k���漰����ɢ�����ĸ�������k����ɢ���㣬k���ӽ�����ɢ����


	void Messagedis();          //��Ϣ�ַ�����
	Message Dmess(int k);                     //����ֵ�ǣ��ַ���Ϣʱ����k�ĸ��Ŵ�����k����Ϣpotential
    vector<Dismes> Dgenbycv(vector<int> & vec);  //����ֵ�ǣ��������ܶȽ�����ɢ���������������׵㼯vec��������ɢ��Ϣ
    //���������漰���������ʵ�ʱ�򽫻�����


    Conpost Comconp(int k, Message & mes);//����������k�ĺ���ֲ�
    Dispost Comdisp(int k, Message & mes);//������ɢ��k�ĺ���ֲ�
	void Composterior();        //����ÿһ����ĺ���ֲ�,��ʱʹ��VE��VR��Ԫû��������ΪVE��Ԫ�����޷����⣬���ԣ��������ɢ��������ʹ��VE��Ԫ
    double Multip(vector<int> & vec, vector<int> & vecs, vector<Dismes> & dmes); //����ֵ�ǣ����漰�Ӻ͵���ɢ������Ϣdmes�У�
    //����dmes�����е���ɢ�������е���ɢ����ȫ��vec�ڹ̶�״̬vecs�£�����ɢ������ȡֵ���г˻�


    double Multiplep(vector<int> inpos, vector<Dismes> & dmes);
    vector<int> Inpos(vector<vector<int> > & allsta, vector<Dismes> & dmes);

    Message Elimdisone(int k, Message & mess); //����������ɢ����k�����µ���Ϣ
    vector<int> Allv(Message & mes);  //������Ϣmes���漰�����е㡣
    vector<int> Setminus(vector<int> & allv, vector<int> & parv); //����allvȥ��parv��ĵ㼯��

    vector<vector<int> > Disstru(vector<Dismes> & dmes); //����dmes�з�֤����ɢ����������ṹ
    vector<int> Elimordering(vector<int> & comdisv, vector<Dismes> & dmes);//comdisv�е���ȫ��������dmes�еı���������Ԫ��
    int Maxpoint(vector<int> & comdisv, vector<int> & elimv);//����elimv��Ȩ�����ĵ㣬
    //���ж��Ȩ�����ĵ���comdisv���е�Ȩ�����ʱ���ȷ���comdisv�е�,comdisv��elimv���Ӽ�
    int Changew(int k, int l, vector<int> & elimv);   //�����Ƿ����Ȩ��


	void Totalmess();   //����������ϵ�����Ϣ������ԭʼ��Ϣ���ɸ��ŷַ���������Ϣ���Լ������Ŵ�������Ϣ
    void Output(vector<Dispost> & dis, vector<Conpost> & con);  //���ı���ʽ������յĽ��
    void Outputstru(int Ver, double cpro, int evidnum, int net_id, std::mt19937& gen);
	void OutputNet(int Ver, double cpro, int evidnum, int net_id, std::mt19937& gen);
	void OutputData(int Ver, double cpro, int evidnum, int n_sample, int net_id, std::mt19937& gen);
	~Propagation();
private:
	vector<Message> omess;   //omess[k]��ʾ��k+1�ϵ�ԭʼ��Ϣpotential
	vector<Message> smess1;  //smess1[k]��ʾ��k+1�ϴ������ĸ��ŵ���Ϣpotential,smess1[0]��ʾ��1�ϵ�ԭʼ��Ϣpotential
	vector<Message> smess2;  //smess2[k]��ʾ����k+1�ĸ��Ŵ�����k+1����Ϣpotential,smess2[0]��ʾ��1�ϵ�ԭʼ��Ϣpotential
    vector<int> meslabel;    //meslabel[k]=0��ʾ��k+1��û�д����ĸ��Ž�����Ϣ��meslabel[k]=1��ʾ��k+1�Ѿ������ĸ��Ž�����Ϣ
    vector<Message> totalmess; //totalmess[k]��ʾ��k+1��������Ϣ������ԭʼ��Ϣ���Ӹ��Ŵ�������Ϣ�Լ������Ŵ�������Ϣ

    vector<int> verfapos;            //verfapos[k]��ʾһ���źţ����ǰ�������k+1���丸����һ���ŵ��ź�
	vector<vector<int> > cliquefun;  //cliquefun[k]��ʾһ�����㼯�������ź�Ϊk+1�����б����ڵ���������,�����ܶȵ��׵㼯

    //����7������������exchange�������ı�
	vector<vector<int> > tverChi;     //��ʱ�Ľ�㼰���ӽ�㼯��
	vector<vector<int> > tverFa;      //��ʱ�ļ��弯��
    vector<vector<int> > tverdisFa;   //��ʱ��V������ļ�������ɢ������
    vector<vector<int> > tverconFa;   //��ʱ��V������ļ���������������
	vector<vector<int> > tstate;   //��ʱ����ɢ״̬�б�
	vector<vector<double> > tprob; //��ʱ�������ֲ�
	vector<int> tdistrnum;         //��ʱ�ĸ��������漰�ķֲ�����

	vector<vector<int> > estate;   //������ɢ֤�ݱ���ȡֵ�����ɢ״̬�б�
	vector<vector<double> > eprob; //������ɢ֤�ݱ���ȡֵ��������ֲ�
	vector<int> edistrnum;         //������ɢ֤�ݱ���ȡֵ��ĸ��������漰�ķֲ���������3������Propagation()����ʱȷ������

    vector<int> desvisit; //��������ȱ���������ʱ����¼�Ƿ��Ѿ�������������㣬1��ʾ�Ѿ�������0��ʾû�б���
    vector<int> dreach; //��¼���Ƿ��Ѿ�dreach,1��ʾdreach,0��ʾû��
    vector<vector<int> > Favisit;//�Ӹ����Ƿ���ʹ���1��ʾ���ʹ���0��ʾû��
    vector<vector<int> > Chivisit;//���ӱ��Ƿ���ʹ���1��ʾ���ʹ���0��ʾû��

    vector<int> evid;        //֤�ݱ�����
    vector<int> poster;      //��֤�ݱ���������Ҫ��������ÿ����ĺ������

     vector<int> weight;           //����ĵ�ǰȨ��
	 vector<int> tweight;          //�������ʱȨ��
	 vector<int> label;            //������(1)��(0)�Ѿ������
	 vector<int> visit;            //������(1)��(0)�����ʹ�
	 vector<vector<int> > stru;    //�����ṹ

};


vector<int> Propagation::Dreach(vector<int> & S, vector<int> & vec){
//       cout<<"����Dreach: "<<endl;
       int i,j,k;
       vector<int> dreachve;
       for(i=0;i<V;i++){dreach[i]=0;}
       Favisit=tverFa; //ʹ��Dreachʱ���Ƚ�Favisit��Chivisit����
       Chivisit=tverChi;
       for(i=0;i<Favisit.size();i++){for(j=0;j<Favisit[i].size();j++){Favisit[i][j]=0;}}
       for(i=0;i<Chivisit.size();i++){for(j=0;j<Chivisit[i].size();j++){Chivisit[i][j]=0;}}
       for(j=0;j<S.size();j++){
          k=S[j];
          dreach[k-1]=1;//��ΪS�еĵ�������d�ɴ���һ��ƽ�������
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
       //cout<<"dreach�Ķ���"<<endl;
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
    vector<int> vect; //vect����l��l�ĺ�����Ĳ���
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
		iter=min_element(vec2.begin(),vec2.end());       //ȡ����������С�ı��
		k=morderver[*iter-1];                            //��С��Ŷ�Ӧ�Ķ���
	    verfapos.push_back(mverpos[k-1]);
	}

    //����ķ��䷽ʽ��Ҫ��д
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
   //cout<<"��ʼ��Ϣ�ռ��� "<<endl;
   int i;
   for(i=0;i<adjcliques[0].size();i++){
	   if(adjcliques[0][i]!=cliquefa[0]){
		   //cout<<"��1������ "<<adjcliques[0][i]<<" ������1����Ϣ"<<endl;
		   Cmess(adjcliques[0][i]);
	   }
   }
   //cout<<"��Ϣ�ռ���� "<<endl;
}


void Propagation::Localstru(Message & mes){
	//cout<<"����Localstru() "<<endl;
	int i,j,m;
    vector<int> vec1;
	vector<vector<int> > vec;         //vec����Ϣ�м���ļ���
	for(i=0;i<mes.cmes.size();i++){
		vec.push_back(mes.cmes[i].cverFa);
	}
	for(i=0;i<mes.dmes.size();i++){
        vec.push_back(mes.dmes[i].dver);
    }

	for(i=0;i<V;i++){tverFa[i].clear();tverFa[i].push_back(i+1);}
	for(i=0;i<vec.size();i++){ //ȷ�����弯
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

	vector<Conmes> conmes=mes.cmes;  //����conmes��������Ϣ��������
	vector<Dismes> dismes=mes.dmes;   //����dismes����ɢ��Ϣ��������
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
	//cout<<"����Localstru()��� "<<endl<<endl;
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


Message Propagation::Cmess(int k){ //����ֵ�ǣ��ռ���Ϣʱ����k�����ŵ���Ϣpotential
	int i,j;
	Message mes=omess[k-1];       //mes��ʾ��k�ϵĵ�ǰ��Ϣ
	vector<int> vec=adjcliques[k-1];
	Message tmes,cmes;
	if(vec.size()==1){
	   cmes=Projmes(mes,k,cliquefa[k-1]);
	}
    else{
		for(i=0;i<vec.size();i++){
		//	cout<<"�� "<<k<<" ������ "<<vec[i]<<" "<<endl;
			if(vec[i]!=cliquefa[k-1]){
		//		cout<<"�� "<<k<<" ������"<<vec[i]<<" ������ "<<k<<" ����Ϣ "<<endl;
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
	smess1[k-1]=cmes;         //��¼����Ϣ�ռ�ʱ����k�������ĸ��ŵ���Ϣ
    return cmes;
}

Message Propagation::Dmess(int k){        //kӦ�ò�Ϊ1,����������1�ĸ���0������1����Ϣ
	//cout<<"���ú���Dmess()"<<endl;
	int i,j,m;
	m=cliquefa[k-1];             //��k�ĸ�����m
	Message mes=omess[m-1];      //��m�ϴ��е�ԭʼ��Ϣ

	vector<int> vec=adjcliques[m-1];  //��m��������
	Message tmes,cmes;
//	cout<<"��Ϣ�ַ����̿�ʼ��"<<endl;
//	cout<<"�� "<<m<<" ���� "<<k<<" �ַ���Ϣ"<<endl;
	for(i=0;i<vec.size();i++){
		if((vec[i]!=k)&&(vec[i]!=0)){	//����m��������vec[i]�Ȳ�����0�����������1�ĸ��ţ��ֲ�����k
//			cout<<"��Ϣ�ַ�ʱ"<<"�� "<<m<<" ������ "<<vec[i]<<" ��������Ϣ"<<endl;
			if(cliquefa[vec[i]-1]==m){
//				cout<<"��Ϣ�ַ�ʱ"<<"�� "<<m<<" ������ "<<vec[i]<<" ��������Ϣ�Ѿ�����Ϣ�ռ�ʱ����¼"<<endl;
			    tmes=smess1[vec[i]-1];

			}
			else {
//				cout<<"��Ϊ�� "<<vec[i]<<" ���� "<<m<<" �ĸ���,��һ��������Ϣ�ַ�"<<endl;
				if(meslabel[m-1]==0){
//					cout<<"�� "<<m<<" δ���ܸ��Ŵ�������Ϣ"<<endl;
					tmes=Dmess(m);
				}
				else{
//					cout<<"�� "<<m<<" ��֮ǰ��Ϣ�ַ������У��Ѿ������˸��� "<<vec[i]<<" ��������Ϣ"<<endl;
					tmes=smess2[m-1];
				}

     /*           cout<<"������Ϣ2: "<<endl;
				for(j=0;j<tmes.cmes.size();j++){
					for(int q=0;q<tmes.cmes[j].cverFa.size();q++){
					    cout<<tmes.cmes[j].cverFa[q]<<" ";
					}cout<<endl;
				}cout<<endl;
				cout<<"��ɢ��Ϣ2: "<<endl;
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
    smess2[k-1]=cmes;    //��¼�ַ���Ϣʱ����k�����ĸ��Ž��ܵ���Ϣ
	meslabel[k-1]=1;
	//cout<<"���ú���Dmess( "<<k<<" )���"<<endl;
    return cmes;
}


Message Propagation::Projmes(Message & mes, int k, int l){
	//cout<<"���ú���Projmes()��"<<endl;
	//��A��ʾ��k�ϵĵ�ǰ��Ϣ�㼯��һ�����������ܶȵ��׵㣬һ�������������ʵ��׵�
	vector<int> A,R,RE,RN,RNC,RND,RNDE;
	vector<int>::iterator it;
    int i,j,m,num;
	vector<int> S=Sep(k, l);
	A=Allv(mes);R=A;
	for(i=0;i<S.size();i++){it=find(R.begin(),R.end(),S[i]);if(it!=R.end()){R.erase(it);}}
	for(i=0;i<R.size();i++){if(elabel[R[i]-1]==0) RN.push_back(R[i]);else RE.push_back(R[i]);}
	Localstru(mes);
	Message remess=Elimbarren(S,RN,mes),mess,mess1;//remess��ʾ��ȥbarren���������µ���Ϣ
	Localstru(remess);
	vector<int> relsvec=Dreach(S,RE);//S�ܹ�d�ɴ�ĵ㼯
    mess=Exemess(relsvec,remess);
    Localstru(mess);


    RN.clear();A=Allv(mess);R=A;
   	for(i=0;i<S.size();i++){it=find(R.begin(),R.end(),S[i]);if(it!=R.end()){R.erase(it);}}
	for(i=0;i<R.size();i++){if(elabel[R[i]-1]==0) RN.push_back(R[i]);}
    for(i=0;i<RN.size();i++){if(dcver[RN[i]-1]==-1) RND.push_back(RN[i]);else RNC.push_back(RN[i]);}

    for(i=RNC.size()-1;i>=0;i--){num=Domnum(RNC[i]);if(num<14) mess=Exchange(RNC[i],mess);}
	//RNC�Ѿ�����С�������򣬴���Ŵ�ĵ㿪ʼexchange��ȥ

	//�в��еķ�֤����ɢ��������û������������Ϊ�ӽ�����ɢ�������Ա�exchange��
	for(i=0;i<RND.size();i++){
            m=0;
            for(j=1;j<tverChi[RND[i]-1].size();j++){
                      if(dcver[tverChi[RND[i]-1][j]-1]==-2) {m++; break;}
            }
            if(m==0) RNDE.push_back(RND[i]);
    }
	for(i=RNDE.size()-1;i>=0;i--){num=Domnum(RNDE[i]);if(num<14) mess=Exchange(RNDE[i],mess);}
    //cout<<"���ú���Projmes()���"<<endl;
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
	vector<int> inpos;         //����ʮ���ƺͶ���Ƶ�ת����ϵ
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
        vector<int> remk=Setminus(ved,vec); //remk�����ֻ��һ��Ԫ��k
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




vector<Conmes> Propagation::Exchange(int k, Conmes & conk, int l, Conmes & conl){ //�ں�������״̬�У�����������k������tverChi����������С����������l����exchange
    //cout<<"���ú���Exchange("<<k<<" , "<<l<<")��������"<<endl;
	int i,j;
	Conmes Conm;
	vector<Conmes>  Vecconm;
	vector<int> veA,veB;
	vector<double> fAk,fBl,fk,fl;
	vector<int> vecd1,vecd2,vecc1,vecc2;
	vector<vector<int> > st1,st2,st3;     //exchange��k��l�й�ͬ����ɢ״̬st3
	vector<vector<double> > pr1,pr2,prk,prl;
	vector<int>::iterator it;
	st1=conk.cverFasta;
	st2=conl.cverFasta;
	pr1=conk.cverFaf;
	pr2=conl.cverFaf;
	vector<int> D(tverdisFa[k-1].size()+tverdisFa[l-1].size());
	vecd1=tverdisFa[k-1];      //��ʱ�̶�סk��l����ɢ����
	vecd2=tverdisFa[l-1];
	merge(vecd1.begin(),vecd1.end(),vecd2.begin(),vecd2.end(),D.begin()); //��k��l����ɢ���弯������
    sort(D.begin(),D.end());
	it=unique(D.begin(),D.end());
	D.erase(it,D.end());                 //ȥ��D���ظ��ĵ�
	st3=EDisstate(D);
    vector<int> C(tverconFa[k-1].size()+tverconFa[l-1].size());
    vecc1=tverconFa[k-1];       //��ʱ�̶�סk��l����������
    vecc2=tverconFa[l-1];
	sort(vecc1.begin(),vecc1.end());
	sort(vecc2.begin(),vecc2.end());
	merge(vecc1.begin(),vecc1.end(),vecc2.begin(),vecc2.end(),C.begin()); //��k��l���������弯������
	sort(C.begin(),C.end());
    it=unique(C.begin(),C.end());
	C.erase(it,C.end());                //ȥ��C���ظ��ĵ�
	it=find(C.begin(),C.end(),k);
	if(it!=C.end()) C.erase(it);        //Cʵ������exchange��k����������㼯
	vector<int> CC=C;//CC��exchange��k��l��ͬ����������㼯������C��ȥ����l�ɵ�CC
	it=find(CC.begin(),CC.end(),l);
	if(it!=CC.end()) CC.erase(it);
    int pos=Pos(k,tverconFa[l-1]);  //��¼��k��l���������弯�е�λ��
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
        fl.push_back(fBl[0]+cof*fAk[0]); //exchange��l�Ļع鳣��
		//exchange��k�Ļع鳣��
		fk.push_back((fAk[0]*fBl[fBl.size()-1]-fBl[0]*cof*fAk[fAk.size()-1])/(fBl[fBl.size()-1]+cof*cof*fAk[fAk.size()-1]));
		for(j=0;j<CC.size();j++){
			if(find(vecc1.begin(),vecc1.end(),CC[j])!=vecc1.end()) a=fAk[Pos(CC[j],tverconFa[k-1])-1];//a=fAk[Tpos(CC[j],k)-1];
			else a=0;//
			if(find(vecc2.begin(),vecc2.end(),CC[j])!=vecc2.end()) b=fBl[Pos(CC[j],tverconFa[l-1])-1];//b=fBl[Tpos(CC[j],l)-1];
			else b=0; //
            fl.push_back(b+cof*a); //exchange��l����CC=C\{l}�ع�ϵ��
			//exchange��k����CC=C\{l}�ع�ϵ��
			fk.push_back((a*fBl[fBl.size()-1]-b*cof*fAk[fAk.size()-1])/(fBl[fBl.size()-1]+cof*cof*fAk[fAk.size()-1]));
		}

		fl.push_back(fBl[fBl.size()-1]+cof*cof*fAk[fAk.size()-1]); //exchange��l�ķ���
		prl.push_back(fl);
		fk.push_back(cof*fAk[fAk.size()-1]/(fBl[fBl.size()-1]+cof*cof*fAk[fAk.size()-1])); //exchange��k����l�Ļع�ϵ��
		fk.push_back(fAk[fAk.size()-1]*fBl[fBl.size()-1]/(fBl[fBl.size()-1]+cof*cof*fAk[fAk.size()-1])); //exchange��k�ķ���
		prk.push_back(fk);
	}

    //�Զ�̬�ĸ��������漰�ķֲ���������̬��ɢ״̬�б���̬�����ֲ������иı�
	tdistrnum[k-1]=st3.size();
	tdistrnum[l-1]=st3.size();
	tstate[k-1].clear();
	tstate[l-1].clear();
	for(i=0;i<st3.size();i++){for(j=0;j<st3[i].size();j++){tstate[k-1].push_back(st3[i][j]);tstate[l-1].push_back(st3[i][j]);}}
	tprob[k-1].clear();
	for(i=0;i<prk.size();i++){for(j=0;j<prk[i].size();j++){tprob[k-1].push_back(prk[i][j]);}}
	tprob[l-1].clear();
	for(i=0;i<prl.size();i++){for(j=0;j<prl[i].size();j++){tprob[l-1].push_back(prl[i][j]);}}

	tverdisFa[k-1]=D;                //���exchange���tverdisFa[k-1],tverdisFa[k-1]�����˸ı�
	tverconFa[k-1].clear();
	tverconFa[k-1].push_back(k);
	for(i=0;i<C.size();i++){tverconFa[k-1].push_back(C[i]);} //���exchange���tverconFa[k-1]��tverconFa[k-1]�����˸ı�

	vector<int> Uk(D.size()+C.size());
	merge(D.begin(),D.end(),C.begin(),C.end(),Uk.begin());
	sort(Uk.begin(),Uk.end());
	tverFa[k-1].clear();
    tverFa[k-1].push_back(k);
	for(i=0;i<Uk.size();i++){tverFa[k-1].push_back(Uk[i]);}//���exchange���tverFa[k-1],tverFa[k-1]�����˸ı�
    it=find(tverChi[k-1].begin(),tverChi[k-1].end(),l);
	tverChi[k-1].erase(it);    //��k���ӽ�㼯��ȥ��l,���exchange���tverChi[k-1]��tverChi[k-1]�����˸ı�
	it=tverChi[l-1].begin();   //��l���ӽ�㼯������k,���exchange���tverChi[l-1]��tverChi[l-1]�����˸ı�
	it++;
	while(it!=tverChi[l-1].end()){if(*it>k){tverChi[l-1].insert(it,k);break;}it++;}
	if(it==tverChi[l-1].end()) tverChi[l-1].push_back(k);
	tverdisFa[l-1]=D;              //���exchange���tverdisFa[l-1]��tverdisFa[l-1]�����˸ı�

	tverconFa[l-1].clear();
	tverconFa[l-1].push_back(l);
	for(i=0;i<CC.size();i++){tverconFa[l-1].push_back(CC[i]);}//���exchange���tverconFa[l-1]��tverconFa[l-1]�����˸ı�
	vector<int> Ul(D.size()+CC.size());
	merge(D.begin(),D.end(),CC.begin(),CC.end(),Ul.begin());
	sort(Ul.begin(),Ul.end());
	tverFa[l-1].clear();
    tverFa[l-1].push_back(l);
	for(i=0;i<Ul.size();i++){tverFa[l-1].push_back(Ul[i]);}//���exchange���tverFa[l-1]��tverFa[l-1]�����˸ı�
	int v;
	for(i=1;i<tverFa[l-1].size();i++){  //exchange���l�ĸ������ӽ�㼯�����˸ı�
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

        vector<Dismes> Vecdism;  //Vecdism[0]��������exchange��k����ɢ��Ϣ��Vecdism[1]����exchange��l����ɢ��Ϣ
        Dismes dism;
        vector<int>::iterator it;
        vector<vector<int> > sk,sl, stk,stl; //stk��stl�ֱ��ʾexchange���k�͵�l�漰����ɢ״̬
        vector<double> pk,pl, prk,prl;//prk��prl�ֱ��ʾexchange���k�͵�l�漰�Ķ�Ӧ״̬�µĸ���
        vector<int> Fa,lFa,kFa;     //Fa��exchange��k��l��ͬ�ĸ�����������ɢ������kFa��exchange��k�ļ��弯������l
        vector<int> posk,posl,veA,veB,veC;
        double p;
        int posA;

        for(i=1;i<tverFa[k-1].size();i++){Fa.push_back(tverFa[k-1][i]);}
        for(i=1;i<tverFa[l-1].size();i++){if(tverFa[l-1][i]!=k) Fa.push_back(tverFa[l-1][i]);}
        sort(Fa.begin(),Fa.end());
        it=unique(Fa.begin(),Fa.end());
	    Fa.erase(it,Fa.end()); //ȥ��Fa�ж���ĵ�
	    lFa.push_back(l);
	    for(i=0;i<Fa.size();i++){lFa.push_back(Fa[i]);}//exchange��l�ļ��弯
	    Fa.push_back(l);  //��Fa����ӵ�l
	    sort(Fa.begin(),Fa.end());
	    kFa.push_back(k);
	    for(i=0;i<Fa.size();i++){kFa.push_back(Fa[i]);} //exchange��k�ļ��弯
	    stk=EDisstate(kFa);//exchange��k�ļ��弯��Ӧ��״̬
	    stl=EDisstate(lFa);//exchange��l�ļ��弯��Ӧ��״̬

	    sk=disk.dversta;//exchangeǰk�ļ��弯��Ӧ��״̬
	    sl=disl.dversta;//exchangeǰl�ļ��弯��Ӧ��״̬
	    pk=disk.dverp;//exchangeǰk�ĸ��ʷֲ�
	    pl=disl.dverp;//exchangeǰl�ĸ��ʷֲ�
        posA=Pos(k,disl.dver);
	    for(i=0;i<stl.size();i++){
            p=0;
            posl=Commpos(lFa,stl[i],disl.dver);  //posl��posk��size()һ��,���ǵ�size()Ӧ�Ǳ���k��ȡֵ����
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
	    tverChi[k-1].erase(it);    //��k���ӽ�㼯��ȥ��l,���exchange���tverChi[k-1]��tverChi[k-1]�����˸ı�
	    it=tverChi[l-1].begin();   //��l���ӽ�㼯������k,���exchange���tverChi[l-1]��tverChi[l-1]�����˸ı�
	    it++;
	    while(it!=tverChi[l-1].end()){if(*it>k){tverChi[l-1].insert(it,k);break;}it++;}
	    if(it==tverChi[l-1].end()) tverChi[l-1].push_back(k);

	    int v;
	    for(i=1;i<tverFa[l-1].size();i++){  //exchange���l�ĸ������ӽ�㼯�����˸ı�
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
	//cout<<"��ʼ��Ϣ�ַ��� "<<endl;
	int i;
	tverFa=verFa;
	tverChi=verChi;
	tverdisFa=verdisFa;
	tverconFa=verconFa;

	tstate=estate;
    tprob=eprob;
	tdistrnum=edistrnum;
	for(i=2;i<cliquenum+1;i++){
//		cout<<"�� "<<i<<" �����ĸ��Ž�����Ϣ"<<endl;
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
			//�����Ե�������ָʾ�����У���data.csv�е�����һ��ʱ���Զ���ȡ���������趨��
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
	
	// �����ͷ
	for (size_t i = 0; i < verFa.size(); ++i) {
		outfile << "V" << (i + 1) << (i < verFa.size() - 1 ? "," : "\n");
	}
	
	// ��������
	for (size_t n = 0; n < n_sample; ++n) {
		//cout << "��" << n + 1 << "������" << endl;
		std::vector<double> nodesValue;
		for (size_t i = 0; i < verFa.size(); ++i) {
			const auto& node_with_parents = verFa[i];
			vector<int> parentsIndex(node_with_parents.begin() + 1, node_with_parents.end());
			vector<int> contParentsIndex;
			vector<int> parentsState;

			for (size_t idx : parentsIndex) {
				// ���ڵ�û�и��ڵ������
				if (dcver[idx - 1] == -1) { // ֻ������ɢ���ڵ�
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
			// ��ɢ����
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
			// ��������
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

				// ֻ��Ҫ�����������ڵ�
				for (int p = 0; p < contParentsIndex.size(); p++) {
					mean = mean + realParas[p + 1] * nodesValue[parentsIndex[p] - 1]; 
				}
				nodeValue = RandGen::normValue(mean, sqrt(sigma2), gen);
			}
			nodesValue.push_back(nodeValue);
		}

		// �������
		for (size_t i = 0; i < nodesValue.size(); ++i) {
			outfile << nodesValue[i] << (i < nodesValue.size() - 1 ? "," : "\n");
		}
	}
}


Message Propagation::Exchange(int k, Message & mes){
	//cout<<"����exchange������ "<<k<<" "<<endl;
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
	      vecm=Exchange(k,conk,l,conl);  //kӦ��������С��l��exchange
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
    for(i=1;i<tverFa[k-1].size();i++){    //��k����exchange�󣬳�Ϊ�����ĵ㣬�����ڽ�һ��ʹ��k֮��ĵ�ʱ
       l=tverFa[k-1][i];                  //��Щ��ĸ�������ӽ���ж����Ậ��k����k�����¼����ȫexchange��k�ĸ����
	   it=find(tverChi[l-1].begin(),tverChi[l-1].end(),k);
	   if(it!=tverChi[l-1].end()) tverChi[l-1].erase(it);
	}
	for(i=0;i<vecmes.size();i++){if(vecmes[i].cverFa[0]!=k) mess.cmes.push_back(vecmes[i]);}
	for(i=0;i<vedmes.size();i++){if(vedmes[i].dver[0]!=k) mess.dmes.push_back(vedmes[i]);}

    return mess;
}


void Propagation::Totalmess(){
    //cout<<"����Totalmess()"<<endl;
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
                        outfile<<"��"<<dis[j].disver<<endl;
                        for(k=0;k<dis[j].state.size();k++){outfile<<dis[j].state[k]<<" ";}outfile<<endl;
                        for(k=0;k<dis[j].dispr.size();k++){outfile<<dis[j].dispr[k]<<" ";}outfile<<endl;
                        break;
                     }
                 }
             }
             else{
                  for(j=0;j<con.size();j++){
                      if(poster[i]==con[j].conver){
                        outfile<<"��"<<con[j].conver<<endl;
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
    //cout<<"����Composterior"<<endl;
    int i,j,k,m;
    vector<int> cliqvex,allvex,rvex,RN,RE; //���з���ĵ㼯,�������е�ļ���,allvex/cliqvex
    vector<int> econ1,econ2;//rvex������������������cliqvex����������������
    vector<int> edis1,edis2;//rvex����������ɢ������cliqvex����������ɢ����
    vector<int> never,partall,rem;
    vector<int> elimordering;
    vector<int> relvec;
    vector<vector<int> > disstru;
    Message mess,mess1;
    vector<Conmes> conmes;
	vector<Dismes> dismes,nedism;
	vector<Message> seqmes;
    Dispost dp;  //��ɢ����
	Conpost cp;  //��������
	vector<Dispost> Dpost;//��ɢ��ĺ�����Ϣ
	vector<Conpost> Cpost;//������ĺ�����Ϣ
	Totalmess();
    for(i=0;i<cliquefun.size();i++){
       //cout<<"�� "<<i+1<<" �ϵ�����"<<endl;
       cliqvex=cliquefun[i];
	   if(cliqvex.size()>=1){
             mess=totalmess[i];
             allvex=Allv(mess);
             rvex=Setminus(allvex,cliqvex);
             RN.clear();RE.clear();
             for(j=0;j<rvex.size();j++){if(elabel[rvex[j]-1]==0) RN.push_back(rvex[j]);else RE.push_back(rvex[j]);}
	         Localstru(mess);  //�Խṹ��������
             mess=Elimbarren(cliqvex,RN,mess);
             Localstru(mess);
             relvec=Dreach(cliqvex,RE);//��cliqvex�ܹ�dreach�ı�����
             mess=Exemess(relvec,mess);//�õ�����ر�����ص���Ϣ
             Localstru(mess);       //�Խṹ��������
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
    //Output(Dpost,Cpost);  ������еĺ���ֲ�
}


vector<Dismes> Propagation::Dgenbycv(vector<int> & vec){//����ֵ�ǣ��������ܶȽ�����ɢ���������������׵㼯vec��������ɢ��Ϣ
	//cout<<"���ú���Dgenbycv()"<<endl;
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
	//cout<<"���ú���Dgenbycv()���"<<endl;
	return Xcdmes;
}


Conpost Propagation::Comconp(int k, Message & mes){
	//cout<<"����Comconp("<<k<<") "<<endl;
	int i,j;
	Conpost cp;
    vector<int> disfa,confa,allcon,alldis,sta,pos,never;
	double sum=0;
	vector<double> mp,mpp,fp;
	vector<vector<int> > disfasta,alldissta,trucall;
	vector<vector<double> > conprob;
	vector<int> A,R,RN,RE;//A���漰��ȫ��㣬R��Aȥ��k��ĵ㼯��RN��R�з�֤�ݱ�����,RE��R��֤�ݱ���
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
    Message remess=Elimbarren(veck,RN,mes),mess,mess1;//remess��ʾ��ȥbarren���������µ���Ϣ
    Localstru(remess);
	relkvec=Dreach(veck,RE);//��k�ܹ�dreach�ı�����
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
	nedism=Dgenbycv(never);//����exchange�����е�������֤�ݱ����������ܶȻ��Ϊ��ɢ����



	disfa=tverdisFa[k-1];



    disfasta=EDisstate(disfa);



	conprob=Tprob(k);

	for(i=0;i<nedism.size();i++){dism.push_back(nedism[i]);}
	mess1.dmes.clear();mess1.cmes.clear();
	for(i=0;i<dism.size();i++){mess1.dmes.push_back(dism[i]);}
    vector<int> elimordering=Elimordering(disfa,dism);
    int size=disfa.size();
    for(j=elimordering.size()-1;j>=size;j--){mess1=Elimdisone(elimordering[j],mess1);}//��disfa�����ɢ�������Ӻ͵�



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
     vector<int> elimordering;         //���㼯��Ӧ����Ԫ˳��
     weight.clear();tweight.clear();label.clear();visit.clear();//weight,tweight,label,visit�ж���˳����elimv��һ�µ�
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
                 pos=i+1;                     //��¼Ȩ�ظ���Ķ���
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
    int k,l,m,n;         //����ѭ��
	int z;              //��¼���Ȩ�Ķ���
	int u;
	vector<int>::iterator iter; //���������͵���
	for(k=0;k<V+1;k++){
        //������
		//cout<<"ѭ��λ�� "<<k<<endl;
		if(k==0) {
			z=V+1;
			slabel[V]=1;
			//������
			//cout<<"��һ������Ŷ��� "<<V+1<<endl;
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
				//�����飬������ͼ���ʱȨ�صĸ������
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

int OperonG::Updateweight(int k,int l){    //��ͼ�ж���k����ŵ�ʱ��δ��Ŷ���l��Ȩ���Ƿ�����
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
                 z=i+1;                     //��¼Ȩ�ظ���Ķ���
		}
	}
	slabel[z-1]=1;
    //������
	//cout<<"����ŵĶ��� "<<z<<endl;
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


Dispost Propagation::Comdisp(int k, Message & mes){//������ɢ�����ĺ���ֲ�ʱ��mes��ֻ����ɢ��Ϣmes.dmes
	//cout<<"����Comdisp()"<<endl;
	int i,j,m,pos;
	vector<int> allver,ksta,allcon,never;   //allver��mes�����е���ɢ��,ksta�ǵ�k��״̬
	vector<vector<int> > allversta;
	Dispost dp;
	double sum=0;
	vector<double> mp,mpp,fp;
    vector<Dismes> dism;
    vector<int> dispk;
    Message mess1=mes;
    dispk.push_back(k);
    vector<int> elimordering=Elimordering(dispk,mes.dmes);
    for(j=elimordering.size()-1;j>=1;j--){mess1=Elimdisone(elimordering[j],mess1);}//��k�����ɢ�������Ӻ͵�
    dism=mess1.dmes;
//    cout<<"it is too long"<<endl;
    allver=Allv(mess1); //k��allver��Ψһ�ķ�֤����ɢ����������allversta.size()���Ǳ���k��ȡֵ����
	allversta=EDisstate(allver);
	for(i=0;i<allversta.size();i++){mp.push_back(Multip(allver,allversta[i],dism));} //Multip�൱�ڽ�Inpos��Multiplep�ϲ���һ��
	for(i=0;i<valnum[k-1];i++){ksta.push_back(i+1);mpp.push_back(0);}
	pos=Pos(k,allver);
	for(i=0;i<allversta.size();i++){m=allversta[i][pos-1];mpp[m-1]=mp[i];}
    for(i=0;i<mpp.size();i++){sum=sum+mpp[i];}
    for(i=0;i<mpp.size();i++){fp.push_back(mpp[i]/sum);}
	dp.disver=k;
	dp.state=ksta;
	dp.dispr=fp;
	//cout<<"�� "<<k<<" �ĺ������"<<endl;
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
      Dismes gedism;      //���ڶ�k��ͣ��²�������ɢ��Ϣ
      vector<Dismes> redism,hadism; //redism��¼���µ���ɢ��Ϣ��hadism��¼Ϊ������k��Ҫ�������ɢ��Ϣ
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
	  allver.erase(it); //allver�����������k�����ɢ�����漰�ı���
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
//	cout<<"����Multip"<<endl;
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
	vector<int> inpos;         //����ʮ���ƺͶ���Ƶ�ת����ϵ
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
//	cout<<"���Multip"<<endl;
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

		ve.push_back(vec); //��Ϊ����һ�����뺯������ʹ��tstate[k-1].size()=0ʱ������ve.size()=1
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








int OperonG::Complete(vector<int> vec){   //����㼯�ڵ���ͼ�ļ�СM���ǻ�ͼ���Ƿ���ȫ
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

int OperonG::Mverpos(){  //��¼�µ����ŵ�λ�õĹ�ϵ,�����ص���ͼ�ļ�Сm���ǻ�ͼ���ŵĸ���
	 int k;
	 int z;
	 int l=0,m=0;      //m��󽫱�ʾM�ؿ������
	 msignver.clear();
	 mcliquever.clear();
	 for(k=0;k<V;k++){
		 z=morderver[k];     //z�Ǳ��Ϊk+1�ĵ�
		 //������
		 //cout<<"z "<<z<<" "<<" m "<<m<<" mmadj[z-1] "<<mmadj[z-1].size()<<endl;
         if(l<=mmadj[z-1].size()){
			  m++;
	          mverpos[z-1]=m;             //��¼��z���ڵ��ŵ��ź�
			  mcliquever.push_front(z);   //��¼���ſ�ʼ�ĵ�
			  if(l!=0) msignver.push_front(morderver[k-1]); //��¼ǰһ���ŵı�ʶ��
		 }
		 else mverpos[z-1]=m;
         l=mmadj[z-1].size();
	 }
	 msignver.push_front(morderver[V-1]);
     for(k=0;k<V;k++){
	     mverpos[k]=m+1-mverpos[k];         //��k+1���ڵ��ŵ��ź�
	 }
	 //�������
/*     for(k=0;k<m;k++){
	   cout<<msignver[k]<<" ";
	 }*/
	 return m;
}

void OperonG::Constructtree(){         //��������ͼ�ļ�Сm���ǻ�ͼ���ŵ�junction tree
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
			   vec.push_back(mnumbering[*iter-1]);    //vec��¼��ʱ�ĵ����ڼ��е�ı��
		}
        iter=min_element(vec.begin(),vec.end());       //ȡ����С�ı��
        x=morderver[*iter-1];                          //��С��Ŷ�Ӧ�ĵ�
		cliquefa.push_back(mverpos[x-1]);                //������С��Ŷ�Ӧ�ĵ������ŵ��ź�
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
	int k,l;//��ѭ����
	int z;         //��¼���Ȩ�Ķ���
	mmadj.clear();

	for(k=0;k<V;k++){
	   if(k==0) {
		   for(l=0;l<V;l++){
		      if(dcver[l]==-1) break;
		   }
		   if(l!=V){z=l+1; mlabel[z-1]=1;}
		   else{z=1;mlabel[z-1]=1;}
	  //   cout<<"��һ������ŵĶ��� "<<z<<endl;
	   }  //�������ɢ�㣬��һ����ŵ��λ������һ��ɢ�㣬���û�У�ѡ��һ�����������ȱ��
	   else z=Mcompare();
	   mnumbering[z-1]=V-k;
	   morderver.push_front(z);
	   Mupdate(z);
	}
	//������
	int m;
	vector<int> vec;
	for(k=0;k<V;k++){                       //���������ڼ�
		   vec.clear();
	       vector<int>::iterator iter=mtriadj[k].begin();
		   for(;iter!=mtriadj[k].end();iter++){
			    m=*iter;
			    if(mnumbering[m-1]>=mnumbering[k])
		            vec.push_back(m);
		   }
		   mmadj.push_back(vec);
	}

	//�����飬��������ڼ�
/*	cout<<"�����ڼ�: "<<endl;
	 for(k=0;k<V;k++){
	    vector<int>::iterator it=mmadj[k].begin();
	   for(;it!=mmadj[k].end();it++){
	        cout<<*it<<" ";
	   }
	   cout<<endl;
	}*/
}

void OperonG::Mupdate(int k){   //����ͼ�ж���k����ŵ�ʱ���������������Ȩ��
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
                 j=mweight[i];      //j��¼����δ��Ŷ����е�һ������Ȩ��
                 z=i+1;             //z��¼���ǵ�һ��Ȩ������δ��Ŷ���
		}
	}
	for(i=z-1;i<V;i++){
        k=mweight[i];
		if((k==j)&&(dcver[i]==-1)){
			if(mlabel[i]==0){
			//	   cout<<"����ŵ���ɢ���� "<<i+1<<endl;
                   mlabel[i]=1;
				   return i+1;
			}
		}
	}
//	cout<<"����ŵ��������� "<<z<<endl;
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
			  for(m=l+1;m<ver.size();m++){     //�Ѽ�����ȫ��
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
	// ���������Ŀ¼
	std::wstring base = L"results";
	std::vector<std::wstring> dirs = { L"structure", L"data", L"evidence"};
	CreateSubDirs(base, dirs);

	// ģ������
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
				cout << "������: " << n_vertex << " ������������: " << p_continuous << " ֤�ݱ�����: " << n_evidence << endl;
				// Randomly generate 10 Bayesian networks.
				for (size_t b = 0; b < 100; ++b) {
					std::mt19937 gen(b); // ��bΪ���ӵ�Mersenne Twister������
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
