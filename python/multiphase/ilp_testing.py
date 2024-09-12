import networkx as nx
from phase_assignment_delays import get_paths


def test_case_1():
    G = nx.DiGraph()

    # Create nodes with attributes
    G.add_node(1, attr="DFF")
    G.add_node(2, attr="AND")
    G.add_node(3, attr="MRG")
    G.add_node(4, attr="NOT")
    G.add_node(5, attr="OR")
    G.add_node(6, attr="SPL")

    # Create edges
    G.add_edges_from([(1, 3), (3, 2), (4, 6), (6, 5)])

    # Expected result: two separate paths
    expected_result = {
        frozenset([1, 2, 3]): [(1, 3, 2)],
        frozenset([4, 5, 6]): [(4, 6, 5)],
    }

    return G, expected_result


def test_case_2():
    G = nx.DiGraph()

    # Create nodes with attributes
    G.add_node(1, attr="XOR")
    G.add_node(2, attr="AND")
    G.add_node(3, attr="MRG")
    G.add_node(4, attr="SPL")

    # Create edges forming a single path
    G.add_edges_from([(1, 3), (3, 4), (4, 2)])

    # Expected result: a single path
    expected_result = {
        frozenset([1, 2, 3, 4]): [(1, 3, 4, 2)],
    }

    return G, expected_result


def test_case_3():
    G = nx.DiGraph()

    # Create nodes with attributes
    G.add_node(0, attr="NOT")
    G.add_node(1, attr="DFF")
    G.add_node(2, attr="AND")
    G.add_node(3, attr="SPL")
    G.add_node(4, attr="SPL")
    G.add_node(5, attr="OR")
    G.add_node(6, attr="MRG")

    # Create edges for intersecting paths
    G.add_edges_from([(0, 4), (1, 3), (3, 2), (4, 6), (6, 5), (3, 6)])

    # Expected result: all paths are merged because they intersect at node 3
    expected_result = {
        frozenset([0, 1, 2, 3, 4, 5, 6]): [(1, 3, 2), (0, 4, 6, 5), (1, 3, 6, 5)],
    }

    return G, expected_result


def test_case_5():
    G = nx.DiGraph()

    # Create nodes with attributes
    G.add_node(1, attr="DFF")
    G.add_node(2, attr="OR")
    G.add_node(3, attr="SPL")
    G.add_node(4, attr="MRG")
    G.add_node(5, attr="AND")

    # Create edges for multiple overlapping paths
    G.add_edges_from([(1, 3), (3, 2), (1, 4), (4, 5)])

    # Expected result: all paths are merged due to shared nodes 3 and 4
    expected_result = {
        frozenset([1, 2, 3, 4, 5]): [(1, 3, 2), (1, 4, 5)],
    }

    return G, expected_result


def run_tests():
    test_cases = [test_case_1, test_case_2, test_case_3, test_case_4, test_case_5]

    for i, test in enumerate(test_cases, 1):
        G, expected_result = test()
        result = get_paths(G)

        print(f"Test Case {i}: {'Pass' if result == expected_result else 'Fail'}")
        print(f"Expected: {expected_result}")
        print(f"Result:   {result}")
        print("-" * 40)


if __name__ == "__main__":
    # Run all test cases
    run_tests()
