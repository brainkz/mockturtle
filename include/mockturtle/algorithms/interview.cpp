
// undirected graph with multiple connected components
// 

struct Node {

    std::vector<Node*> neighbors;
}

struct Graph {
    std::vector<Node> nodes;
}



// 1. Take a node 
// 2. DFS - std::stack, std::vector
// 3. Once DFS terminates, go to 1.

int main()
{
    std::unordered_set<Node> nodes;

    std::unordered_set<Node> nodes_copy = nodes;


    std::vector<std::vector<Node>> components;
    
    while (!nodes_copy.empty())
    {
        std::vector<Node> stack;
        stack.push_back(nodes_copy.back());
        nodes_copy.pop_back();

        std::vector<Node> component;
        std::unordered_set<Node> seen;

        while (!stack.empty())
        {
            Node node = stack.back();
            stack.pop_back();

            seen.emplace(node);
            components.push_back(node);
            nodes_copy.remove(node); // log(n) removal inefficiency

            for (auto it = node.neighbors.begin(); it != node.neighbors.end(); ++it)
            {
                // node is unseeen
                if (seen.find(*it) != seen.end())
                {
                    stack.push_back(*it);
                }
            }
        }
        components.push_back(component);
    }
}