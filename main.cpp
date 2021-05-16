#include <iostream>
#include <memory>
#include <vector>
#include <list>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include <set>
#include <fstream>
#include <regex>

// Work of Sinii Viacheslav

using namespace std;

namespace exceptions
{
    class ActionWithNonExistingVertexException : public exception
    {
        [[nodiscard]] const char* what() const noexcept override
        {
            return "ActionWithNonExistingVertexException!";
        }
    };
}

// class with info and method for vertex
template <typename T>
class Vertex
{
public:
    Vertex(unique_ptr<T> in_value, int m_ind)
            :
            value(move(in_value)),
            index_in_matrix(m_ind)
    {}

    void printInfo() const
    {
        cout << "Index: " << index_in_matrix
             << " value: " << *value << endl;
    }

public:
    [[nodiscard]] T& getValue() const {
        return *value;
    }

    [[nodiscard]] int getIndexInMatrix() const {
        return index_in_matrix;
    }

private:
    unique_ptr<T> value;
    int index_in_matrix = -1;
};

// class with info and method for edge
template <typename N, typename T>
class Edge
{
public:
    Edge(shared_ptr<Vertex<T>>& out, shared_ptr<Vertex<T>>& in, N& weight)
            :
            origin(out),
            destination(in),
            weight(weight)
    {}

    void printInfo() const
    {
        cout << "Edge: " << origin->getValue()
             << "-" << destination->getValue()
             << " weight: " << weight << endl;
    }
public:
    [[nodiscard]] N getWeight() const {
        return weight;
    }

    void setWeight(N new_weight) {
        weight = new_weight;
    }

    [[nodiscard]] const shared_ptr<Vertex<T>> &getOrigin() const {
        return origin;
    }

    [[nodiscard]] const shared_ptr<Vertex<T>> &getDestination() const {
        return destination;
    }

private:
    N weight;
    N bandwidth;
    shared_ptr<Vertex<T>> origin;
    shared_ptr<Vertex<T>> destination;
};

// graph interface
template <typename T, typename N>
class GraphADT
{
public:
    virtual shared_ptr<Vertex<T>> addVertex(const T&) = 0;
    virtual shared_ptr<Edge<N, T>> addEdge(shared_ptr<Vertex<T>>, shared_ptr<Vertex<T>>, N) = 0;
    [[nodiscard]] virtual pair<bool, list<shared_ptr<Edge<N, T>>>> edgesFrom(const shared_ptr<Vertex<T>>& vertex) const = 0;
    [[nodiscard]] virtual pair<bool, list<shared_ptr<Edge<N, T>>>> edgesTo(const shared_ptr<Vertex<T>>& vertex) const = 0;
    [[nodiscard]] virtual pair<bool, shared_ptr<Vertex<T>>> findVertex(const T&) const = 0;
    [[nodiscard]] virtual pair<bool, shared_ptr<Edge<N, T>>> findEdge(const T& from, const T& to) const = 0;
    [[nodiscard]] virtual bool hasEdge(const shared_ptr<Vertex<T>>&, const shared_ptr<Vertex<T>>&) const = 0;
};

// my graph class
template <typename T, typename N>
class AdjacencyMatrixGraph : GraphADT<T, N> {
public:
    // graph constructor
    AdjacencyMatrixGraph()
            :
            greatest_occupied_index(-1),
            initial_matrix_size(8),
            matrix(8, vector<shared_ptr<Edge<N, T>>>(8, nullptr))
    {
        for (int i = 0; i < initial_matrix_size; ++i)
        {
            free_vertices.push(i);
        }
    }

    /* Add new vertex */
private:
    // assigns the lowest available vertex to the new vertex if matrix is not full
    inline int simpleAdd()
    {
        // extract the lowest available index
        int index = free_vertices.top();
        free_vertices.pop();

        // check if it is new greatest occupied vertex
        if (index > greatest_occupied_index)
            greatest_occupied_index = index;

        // return assigned index
        return index;
    }

    // if the graph is full - we double its size and then return assigned index
    inline int doubleMatrix()
    {
        // assign the value of index
        int index = matrix.size();
        // as it guaranteed will be the highest - update greatest_occupied_index
        greatest_occupied_index = index;

        // double the size of the matrix
        int newSize = matrix.size() * 2;

        // push all newly available vertices to the set of available vertices
        for (int i = greatest_occupied_index + 1; i < newSize; ++i)
        {
            free_vertices.push(i);
        }

        // double the number of rows
        matrix.resize(newSize);
        // double the number of columns
        for (int i = 0; i < newSize; ++i) {
            matrix[i].resize(newSize);
        }

        // return assigned index
        return index;
    }

    // return assigned index
    inline int newElementIndex()
    {
        // simpleAdd() - if we have free indices.
        // Otherwise, call doubleMatrix()
        return !free_vertices.empty() ? simpleAdd() : doubleMatrix();
    }

public:
    // adds new vertex to the graph
    shared_ptr<Vertex<T>> addVertex(const T& value) override
    {
        // what to do duplicate names of vertices
//        cout << findVertex(value).first << " " << findVertex(value).second << endl;
        if (findVertex(value).first) return findVertex(value).second;

        // obtain index for new vertex
        int index = newElementIndex();

        // create new vertex
        auto new_vertex = make_shared<Vertex<T>>(Vertex<T>(make_unique<T>(value), index));
        // insert new vertex to the list of vertices
        auto it = vertexList.insert({ value, new_vertex });

//        printInfo();

        // return reference to new vertex
        return new_vertex;
    }

    /* End add new vertex */

    /* Add new edge */

private:
    // check if the client tries to add edge with non-existing end vertices
    inline void nonExistingEndVerticesCheck(T fromValue, T toValue)
    {
        if (vertexList.find(fromValue) == vertexList.end()
            || vertexList.find(toValue) == vertexList.end())
            throw exceptions::ActionWithNonExistingVertexException();
    }

    // manage adding new edge to the matrix and edge list.
    // check on adding new edge twice and adding edge with non-existing end vertices
    void manageEdgeAddition(shared_ptr<Edge<N, T>> new_edge)
    {
        // obtain end vertices
        auto fromVertex = new_edge->getOrigin();
        auto toVertex = new_edge->getDestination();

        // obtain values of end vertices
        auto fromValue = fromVertex->getValue();
        auto toValue = toVertex->getValue();

        nonExistingEndVerticesCheck(fromValue, toValue);

        // TODO: adding already existing edge
        // if the client tries to add already existing edge - just update its weight
        auto edge = findEdge(fromValue, toValue);
        if (edge.first)
            edge.second->setWeight(new_edge->getWeight());

        // add the edge to the list of edges
        auto it = edgeList.insert({{ fromVertex, toVertex }, new_edge });

        // fill the according place in the matrix
        matrix[fromVertex->getIndexInMatrix()][toVertex->getIndexInMatrix()] = new_edge;
    }

public:
    // add new edge with bandwidth to the graph
    shared_ptr<Edge<N, T>> addEdge(shared_ptr<Vertex<T>> from, shared_ptr<Vertex<T>> to, N weight) override
    {
        // construct new edge
        auto new_edge = make_shared<Edge<N, T>>(
                Edge<N, T>(
                        from, to, weight
                )
        );

        manageEdgeAddition(new_edge);

//        printInfo();

        // return a reference to the newly created edge object
        return new_edge;
    }

    /* End add new edge */

    // obtain the list of edges outgoing from the specified vertex
    pair<bool, list<shared_ptr<Edge<N, T>>>> edgesFrom(const shared_ptr<Vertex<T>>& vertex) const override
    {
        list<shared_ptr<Edge<N, T>>> edges;
        // if a client tries to find edges from non-existing vertex
        if (vertexList.find(vertex->getValue()) == vertexList.end()) return { false, edges };
//            throw exceptions::ActionWithNonExistingVertexException();


        // obtain vertex's index in the matrix
        int index = vertex->getIndexInMatrix();

        // go through each cell in the vertex’s row in the matrix
        // and collect all non-nullptr values into a list
        for (int i = 0; i < greatest_occupied_index + 1; ++i)
        {
            if (matrix[index][i] != nullptr)
            {
                edges.push_back(matrix[index][i]);
            }
        }

        // if there is no outgoing edges - return nullopt
        if (edges.size() == 0) return { false, edges };
        else return { true, edges }; // otherwise, return the list
    }

    // obtain the list of edges incoming to the specified vertex
    pair<bool, list<shared_ptr<Edge<N, T>>>> edgesTo(const shared_ptr<Vertex<T>>& vertex) const override
    {
        list<shared_ptr<Edge<N, T>>> edges;

        if (vertexList.find(vertex->getValue()) == vertexList.end()) return { false, edges };
//            throw exceptions::ActionWithNonExistingVertexException();


        // obtain vertex's index in the matrix
        int index = vertex->getIndexInMatrix();

        // go through each cell in the vertex’s column in the matrix
        // and collect all non-nullptr values into a list
        for (int i = 0; i < greatest_occupied_index + 1; ++i)
        {
            if (matrix[i][index] != nullptr)
            {
                edges.push_back(matrix[i][index]);
            }
        }

        // if there is no incoming edges - return nullopt
        if (edges.size() == 0) return { false, edges };
        else return { true, edges }; // otherwise, return the list
    }

    // find a vertex with the specified value
    pair<bool, shared_ptr<Vertex<T>>> findVertex(const T& value) const override
    {
        // try to find a vertex with the specified value in the vertex list
        auto tmp = vertexList.find(value);

        // if the vertex is found - return a reference to it
        if (tmp != vertexList.end())
            return { true, tmp->second };
        else return { false, nullptr }; // if there is no such vertex - return nullopt
    }

    // find an edge with the specified values at end vertices
    pair<bool, shared_ptr<Edge<N, T>>> findEdge(const T& from, const T& to) const override
    {
        // try to obtain vertices with the specified values in the graph
        const auto v1 = findVertex(from);
        const auto v2 = findVertex(to);

        // if the vertices were found
        if (v1.first && v2.first)
        {
            // try to find an edge with the same end vertices in the edge list
            const auto tmp = edgeList.find({ v1.second, v2.second });

            // if the edge is found - return a reference to it
            if (tmp != edgeList.end())
                return { true, tmp->second };
        }

        return { false, nullptr }; // if there is no such edge - return nullopt
    }

    // check if there is an edge between two vertices
    bool hasEdge(const shared_ptr<Vertex<T>>& from, const shared_ptr<Vertex<T>>& to) const override
    {
        // check if there are such vertices in the graph
        if (vertexList.find(from->getValue()) == vertexList.end() || vertexList.find(to->getValue()) == vertexList.end())
            throw exceptions::ActionWithNonExistingVertexException();

        // true if this function could find the edge in the graph
        return findEdge(from->getValue(), to->getValue()).first;
    }


public:
    // print info about contents of a graph
    void printInfo() const
    {
        cout << "Matrix (size is " << matrix.size() << "): " << endl;
        for (int i = 0; i < greatest_occupied_index + 1; ++i)
        {
            for (int j = 0; j < greatest_occupied_index + 1; ++j)
            {
                if (matrix[i][j] == nullptr) cout << 0 << ' ';
                else cout << 1 << ' ';
            }
            cout << endl;
        }

        cout << "Vertex List: " << endl;
        for (const auto& tmp : vertexList)
        {
            tmp.second->printInfo();
        }
        cout << endl;

        cout << "Edge List: " << endl;
        for (const auto& tmp : edgeList)
        {
            tmp.second->printInfo();
        }
        cout << endl;
    }

private:
    // for enabling hashing pairs
    struct pair_hash
    {
        template <class T1, class T2>
        int operator() (pair<T1, T2> const& v) const
        {
            return std::hash<T1>()(v.first) + std::hash<T2>()(v.second);
        }
    };

private:
    unordered_multimap<pair<shared_ptr<Vertex<T>>, shared_ptr<Vertex<T>>>, shared_ptr<Edge<N, T>>, pair_hash> edgeList;
    unordered_multimap<T, shared_ptr<Vertex<T>>> vertexList;
    vector<vector<shared_ptr<Edge<N, T>>>> matrix;

    // contains indices that can be assigned to a new vertex
    priority_queue<int, vector<int>, greater<>> free_vertices;
    int greatest_occupied_index;

    int initial_matrix_size;
};

class FSA : public AdjacencyMatrixGraph<string, string>
{
private:
    class MalformedFileException : public exception
    {
        [[nodiscard]] const char* what() const noexcept override
        {
            return "E0: Input file is malformed";
        }
    };

    class StateNotInSetException : public exception
    {
    public:
        explicit StateNotInSetException(const string& state)
                :
                out_state(string("E3: A state '")
                          + state
                          + string("' is not in the set of states"))
        {}

        [[nodiscard]] const char* what() const noexcept override
        {
            return out_state.c_str();
        }
    private:
        string out_state;
    };

    class DisjointStatesException : public exception
    {
        [[nodiscard]] const char* what() const noexcept override
        {
            return "E2: Some states are disjoint";
        }
    };

    class TransitionNotInSetException : public exception
    {
    public:
        explicit TransitionNotInSetException(const string& trans)
            :
            out_trans(string("E3: A transition '")
                + trans
                + string("' is not represented in the alphabet"))
        {}

        [[nodiscard]] const char* what() const noexcept override
        {
            return out_trans.c_str();
        }
    private:
        string out_trans;
    };

    class NoInitialStateException : public exception
    {
        [[nodiscard]] const char* what() const noexcept override
        {
            return "E4: Initial state is not defined";
        }
    };

    class NondeterministicFSAException : public exception
    {
        [[nodiscard]] const char* what() const noexcept override
        {
            return "E5: FSA is nondeterministic";
        }
    };

public:
    FSA(const string& filename)
        :
        infoFile(filename)
    {
        parseFile();
        constructFsa();
    }

private:
    string readFile()
    {
        string line;
        string result;
        ifstream inputFile(infoFile);

        if (inputFile.is_open())
        {
            while ( getline (inputFile,line) )
            {
                result += line + '\n';
            }
            inputFile.close();
        }

        inputFile.close();

        return result;
    }

    vector<string> splitLines()
    {
        stringstream fsaInfo(readFile());

        vector<string> lines;
        string tmpLine;

        while (getline(fsaInfo, tmpLine, '\n'))
            lines.push_back(tmpLine);

        return lines;
    }

    [[nodiscard]] vector<string> parseLine(const string& line)
    {
        auto setValues = new vector<string>();
        std::regex rgx("([a-zA-Z0-9_>]+)");
        smatch res;

        string::const_iterator searchStart(line.cbegin());
        bool isFirst = true;
        while (regex_search(searchStart, line.cend(), res, rgx))
        {
            if (!isFirst) setValues->push_back(res.str());
            else isFirst = false;

            searchStart = res.suffix().first;
        }

        return *setValues;
    }

    inline bool checkStates(const string& line)
    {
        return line.substr(0, 8) == "states=["
               && line.back() == ']'
               && line[line.length() - 2] != ',';
    }

    inline bool checkAlpha(const string& line)
    {
        return line.substr(0, 7) == "alpha=["
               && line.back() == ']'
               && line[line.length() - 2] != ',';
    }

    inline bool checkInitial(const string& line)
    {
        return line.substr(0, 9) == "initial=["
               && line.back() == ']'
               && line[line.length() - 2] != ',';
    }

    inline bool checkAccepting(const string& line)
    {
        return line.substr(0, 11) == "accepting=["
               && line.back() == ']'
               && line[line.length() - 2] != ',';
    }

    inline bool checkTrans(const string& line)
    {
        return line.substr(0, 7) == "trans=["
               && line.back() == ']'
               && line[line.length() - 2] != ',';
    }

    void parseFile()
    {
        vector<string> lines(splitLines());

        if (!(lines.size() == 5
              && checkStates(lines[0])
              && checkAlpha(lines[1])
              && checkInitial(lines[2])
              && checkAccepting(lines[3])
              && checkTrans(lines[4])))
            throw MalformedFileException();

        states = parseLine(lines[0]);
        for (const auto& tmp : parseLine(lines[1]))
            alphabet.insert(tmp);
        initialInfo = parseLine(lines[2]);
        acceptingInfo = parseLine(lines[3]);
        trans = parseLine(lines[4]);

        if (states.empty()
                || alphabet.empty()
                || trans.empty())
            throw MalformedFileException();
    }

    inline void addVertices()
    {
        for (const auto& tmp : states)
        {
            addVertex(tmp);
        }
    }

    void addEdges()
    {
        string fromVertex,
               toVertex,
               transName;

        for (const auto& tmp : trans)
        {
            stringstream ss(tmp);

            getline(ss, fromVertex, '>');
            getline(ss, transName, '>');
            getline(ss, toVertex);

            auto v1 = findVertex(fromVertex);
            auto v2 = findVertex(toVertex);

            if (!v1.first) throw StateNotInSetException(fromVertex);
            else if (!v2.first) throw StateNotInSetException(toVertex);
            else if (alphabet.find(transName) == alphabet.end()) throw TransitionNotInSetException(transName);
            else
            {
                if (hasEdge(v1.second, v2.second))
                {
                    auto edge = findEdge(fromVertex, toVertex).second;
                    auto newName = edge->getWeight() + "|" + transName;
                    edge->setWeight(newName);
                }
                else addEdge(v1.second, v2.second, transName);
            }
        }
    }

    inline void checkIfDisjoint()
    {
        for (const auto& tmp : states)
        {
            auto vertex = findVertex(tmp);
            auto from = edgesFrom(vertex.second);
            auto to = edgesTo(vertex.second);

            bool toOther = false;
            bool fromOther = false;
            if (from.first)
            {
                for (const auto& tmp1 : from.second)
                {
                    if (tmp1->getDestination()->getValue() != tmp)
                        toOther = true;
                }
            }

            if (to.first)
            {
                for (const auto& tmp1 : to.second)
                {
                    if (tmp1->getOrigin()->getValue() != tmp)
                        fromOther = true;
                }
            }

            if (!((from.first || to.first) && (fromOther || toOther)))
            {
                throw DisjointStatesException();
            }
        }


    }

    inline void checkIfNoInitial()
    {
        if (initialInfo.empty()) throw NoInitialStateException();
    }

    void checkIfNondeterministic()
    {
        string fromVertex, transName;
        unordered_set<string> help;

        for (const auto& tmp : trans)
        {
            int pos = tmp.find('>', tmp.find('>') + 1);

            help.insert(tmp.substr(0, pos));
        }

        if (help.size() != trans.size()) throw NondeterministicFSAException();

    }

    void constructFsa()
    {
        addVertices();
        addEdges();
        checkIfDisjoint();
        checkIfNoInitial();
        checkIfNondeterministic();

        initial = findVertex(*initialInfo.begin()).second;
        for (const auto& tmp : acceptingInfo)
        {
            accepting.push_back(findVertex(tmp).second);
        }

//        printInfo();
//        for (const auto& tmp : trans)
//            cout << tmp << endl;
    }

    void initTable(vector<vector<string>>& table)
    {
//        for (const auto& tmp : states)
//            cout << tmp << endl;

        for (const auto& state : states)
        {
            auto thisVertex = findVertex(state).second;
            auto edges = edgesFrom(thisVertex).second;
//            cout << "outgoing edges: " << endl;
//            for (const auto& tmp : edges)
//                tmp->printInfo();
//            cout << "In: ";
//            thisVertex->printInfo();

            unordered_set<int> visited;
            for (const auto& edge : edges)
            {
                auto to = edge->getDestination();
                if (table[thisVertex->getIndexInMatrix()][to->getIndexInMatrix()] == "{}")
                    table[thisVertex->getIndexInMatrix()][to->getIndexInMatrix()] = "";

//                cout << "indices: " << thisVertex->getIndexInMatrix() << " " << to->getIndexInMatrix() << endl;
                table[thisVertex->getIndexInMatrix()][to->getIndexInMatrix()]
                    .append(edge->getWeight())
                    .append("|");
                visited.insert(to->getIndexInMatrix());
            }

            for (const auto& tmp : visited)
            {
                table[thisVertex->getIndexInMatrix()][tmp].pop_back();
            }
        }

        for (int i = 0; i < states.size(); ++i) {
            if (table[i][i] == "{}")
                table[i][i] = "eps";
            else
                table[i][i].append("|eps");
        }

//        printTable(table);
    }

    void constructTable(vector<vector<string>>& table)
    {
        int numberOfStates = states.size();
        for (int k = 0; k < numberOfStates; ++k) {
            vector<vector<string>> previous(table);

            for (int i = 0; i < numberOfStates; ++i)
            {
                for (int j = 0; j < numberOfStates; ++j)
                {
                    table[i][j] =     "(" + previous[i][k] + ")"
                                      + "(" + previous[k][k] + ")*"
                                      + "(" + previous[k][j] + ")"
                                      + "|"
                                      + "(" + previous[i][j] + ")";
                }
            }
//            cout << "table" << endl;
//            printTable(table);
        }
//        printTable(table);
    }

    string makeRegex(const vector<vector<string>>& table)
    {
        int initialIndex = initial->getIndexInMatrix();
        string result;

        if (accepting.empty()) result = "";
        else
        {
            for (const auto& tmp : accepting)
            {
                int finalIndex = tmp->getIndexInMatrix();

                result += table[initialIndex][finalIndex] + "|";
            }

            result.pop_back();
        }

        return result;
    }

    void printTable(const vector<vector<string>>& table)
    {
        for (int i = 0; i < table.size(); ++i) {
            for (int j = 0; j < table.size(); ++j) {
                cout << table[i][j] << " ";
            }
            cout << endl;
        }
    }

public:
    string toRegex()
    {
        vector<vector<string>> table(
                states.size(), vector<string>(
                        states.size(), "{}"));
        initTable(table);
        constructTable(table);

        return makeRegex(table);
    }

private:
    vector<string> states;
    unordered_set<string> alphabet;
    vector<string> initialInfo;
    vector<string> acceptingInfo;
    vector<string> trans;

    shared_ptr<Vertex<string>> initial;
    vector<shared_ptr<Vertex<string>>> accepting;

    string infoFile;
};

int main()
{
    ofstream outFile("output.txt");
    try
    {
        FSA myFsa("input.txt");
//        myFsa.printInfo();
        auto result = myFsa.toRegex();

        if (result.empty()) outFile << "{}";
        else outFile << myFsa.toRegex();
    }
    catch (exception& e)
    {
        outFile << "Error: " << endl;
        outFile << e.what() << endl;
    }

    outFile.close();
    return 0;
}
