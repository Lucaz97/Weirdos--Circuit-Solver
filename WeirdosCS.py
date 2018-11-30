import numpy
import matplotlib.pyplot as plot
import sys

class Edge:
    def __init__(self, component, first, second):
        self.node_one = first
        self.node_two = second
        self.component = component  # The component this Edge belongs to
        self.current_tracker = []
        self.voltage_tracker = []


class Component:
    def __init__(self, type, name, nodes, edges, params):
        self.type = type
        self.identifier = name
        self.node_list = nodes
        self.edge_list = edges
        self.parameter_list = params
        if type == "L" or type == "C":
            self.is_dynamic = True
        else:
            self.is_dynamic = False
        self.is_tracked = False


class Topology:
    def __init__(self, nodes, edges):
        self.node_list = []
        self.edge_list = []
        self.init_clean(nodes, edges)

        self.A_matrix = numpy.zeros((len(self.node_list), len(self.edge_list)))
        self.M_matrix = numpy.zeros((len(self.edge_list), len(self.edge_list)))
        self.N_matrix = numpy.zeros((len(self.edge_list), len(self.edge_list)))
        self.Z_vector = numpy.zeros(len(self.edge_list))
        self.T_matrix = []

    def build_all(self, step):
        self.build_A()
        self.build_M(step)
        self.build_N(step)
        self.build_Z(step)
        self.T_matrix = self.build_T()

    def init_clean(self, nodes, edges):  # Clean nodes, clean edges of nonexistent nodes
        nodes_clean = list(set(nodes))  # Remove duplicates in node vector
        for edge in edges:  # Remove pointers to nonexistent nodes in edges
            for node in nodes_clean:
                if edge.node_one == node:
                    edge.node_one = node
                if edge.node_two == node:
                    edge.node_two = node
        self.node_list = nodes_clean
        self.edge_list = edges

    def build_A(self):
        for edge in self.edge_list:
            self.A_matrix[self.node_list.index(edge.node_one)][self.edge_list.index(edge)] = 1
            self.A_matrix[self.node_list.index(edge.node_two)][self.edge_list.index(edge)] = -1
        # check matrix A consistency
        absM = numpy.absolute(self.A_matrix)
        (l, n) = self.A_matrix.shape
        # check sum of numbers in columns SHOULD BE 0
        rowSum = numpy.sum(self.A_matrix, axis=0)
        for ints in rowSum:
            if ints != 0.0:
                sys.exit("Topology error: Matrix A check failed, at least one edge is wrongly defined. ")
        # check hanging and isolated nodes
        colSum = numpy.sum(absM, axis=1)
        for ints in colSum:
            if ints <= 1:
                sys.exit("Topology error: Matrix A check failed, isolated node detected. ")
        # matrix A is well constructed, proceed to reduce it

        # We're looking for the reference node
        reference_index = -1
        reference_count = 0
        for node in self.node_list:
            if node == "GND":
                reference_index = self.node_list.index(node)
                reference_count += 1
        if reference_count != 1:
            sys.exit("Topology error: couldn't resolve reference node.")

        self.A_matrix = numpy.delete(self.A_matrix, reference_index, axis=0)

    def build_M(self, h):
        for edge in self.edge_list:
            # WARNING: VERY COMPONENT-SPECIFIC
            component = edge.component
            if component.type == "R":  # If component is a Resistor
                self.M_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = 1
            elif component.type == "DCVS" or \
                 component.type == "ACVS" or \
                 component.type == "SWVS" or \
                 component.type == "TWVS":  # If component is a Voltage Source
                self.M_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = 1
            elif component.type == "DCCS" or \
                 component.type == "ACCS" or \
                 component.type == "SWCS" or \
                 component.type == "TWCS":  # If component is a Current Source
                pass
            elif component.type == "L":  # If component is an Inductor
                if self.edge_list.index(edge) == self.edge_list.index(component.edge_list[0]):
                    self.M_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = 1
                else:
                    pass
            elif component.type == "C":  # If component is a Capacitor
                if self.edge_list.index(edge) == self.edge_list.index(component.edge_list[0]):
                    self.M_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = 1
                else:
                    self.M_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = 1
            elif component.type == "DBR":  # If component is Double Bipole with R Matrix
                if component.edge_list.index(edge) == 0:
                    self.M_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = 1
                elif component.edge_list.index(edge) == 1:
                    self.M_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = 1
                else:
                    pass
            elif component.type == "DBG":  # If component is Double Bipole with G Matrix
                if component.edge_list.index(edge) == 0:
                    self.M_matrix[self.edge_list.index(component.edge_list[0])][self.edge_list.index(edge)] = -1 * \
                        float(component.parameter_list[0])
                    self.M_matrix[self.edge_list.index(component.edge_list[1])][self.edge_list.index(edge)] = -1 * \
                        float(component.parameter_list[1])
                elif component.edge_list.index(edge) == 1:
                    self.M_matrix[self.edge_list.index(component.edge_list[0])][self.edge_list.index(edge)] = -1 * \
                        float(component.parameter_list[2])
                    self.M_matrix[self.edge_list.index(component.edge_list[1])][self.edge_list.index(edge)] = -1 * \
                        float(component.parameter_list[3])
                else:
                    pass
            elif component.type == "DBH":  # If component is Double Bipole with H Matrix
                if component.edge_list.index(edge) == 0:
                    self.M_matrix[self.edge_list.index(component.edge_list[0])][self.edge_list.index(edge)] = 1
                    self.M_matrix[self.edge_list.index(component.edge_list[0])][self.edge_list.index( component.edge_list[1])] = -1 *\
                        float(component.parameter_list[1])
                elif component.edge_list.index(edge) == 1:
                    self.M_matrix[self.edge_list.index(component.edge_list[1])][self.edge_list.index(edge)] = -1 * \
                        float(component.parameter_list[3])
                else:
                    pass
            elif component.type == "DBh":  # If component is Double Bipole with H prime Matrix
                if component.edge_list.index(edge) == 0:
                    self.M_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = -1 *\
                        float(component.parameter_list[0])
                elif component.edge_list.index(edge) == 1:
                    self.M_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = 1
                    self.M_matrix[self.edge_list.index(edge)][self.edge_list.index(component.edge_list[0])] = -1 * \
                        float(component.parameter_list[2])
                else:
                    pass
            elif component.type == "DBT":  # If component is Double Bipole with T Matrix
                if component.edge_list.index(edge) == 0:
                    self.M_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = 1
                    self.M_matrix[self.edge_list.index(edge)][self.edge_list.index(component.edge_list[1])] = -1 * \
                        float(component.parameter_list[0])
                elif component.edge_list.index(edge) == 1:
                    self.M_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = -1 * \
                        float(component.parameter_list[1])
                else:
                    pass
            else:
                pass

    def build_N(self, h):
        for edge in self.edge_list:
            # WARNING: VERY COMPONENT-SPECIFIC
            component = edge.component
            if component.type == "R":  # If component is a Resistor
                self.N_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = -1 * float(component.parameter_list[0])
            elif component.type == "DCVS" or \
                 component.type == "ACVS" or \
                 component.type == "SWVS" or \
                 component.type == "TWVS":  # If component is a Voltage Source
                pass
            elif component.type == "DCCS" or \
                 component.type == "ACCS" or \
                 component.type == "SWCS" or \
                 component.type == "TWCS":  # If component is a Current Source
                self.N_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = 1
            elif component.type == "L":  # If component is an Inductor
                if self.edge_list.index(edge) == self.edge_list.index(component.edge_list[0]):
                    self.N_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = \
                        -1 * float(component.parameter_list[0]) / h
                else:
                    self.N_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = 1
            elif component.type == "C":  # If component is a Capacitor
                if self.edge_list.index(edge) == self.edge_list.index(component.edge_list[0]):
                    self.N_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] =\
                        -1 * h / float(component.parameter_list[0])
                else:
                    pass
            elif component.type == "DBR":  # If component is Double Bipole with R Matrix
                if component.edge_list.index(edge) == 0:
                    self.N_matrix[self.edge_list.index(component.edge_list[0])][self.edge_list.index(edge)] = -1 * \
                        float(component.parameter_list[0])
                    self.N_matrix[self.edge_list.index(component.edge_list[1])][self.edge_list.index(edge)] = -1 * \
                        float(component.parameter_list[1])
                elif component.edge_list.index(edge) == 1:
                    self.N_matrix[self.edge_list.index(component.edge_list[0])][self.edge_list.index(edge)] = -1 * \
                        float(component.parameter_list[2])
                    self.N_matrix[self.edge_list.index(component.edge_list[1])][self.edge_list.index(edge)] = -1 * \
                        float(component.parameter_list[3])
                else:
                    pass
            elif component.type == "DBG":  # If component is Double Bipole with G Matrix
                if component.edge_list.index(edge) == 0:
                    self.N_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = 1
                elif component.edge_list.index(edge) == 1:
                    self.N_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = 1
                else:
                    pass
            elif component.type == "DBH":  # If component is Double Bipole with H Matrix
                if component.edge_list.index(edge) == 0:
                    self.N_matrix[self.edge_list.index(component.edge_list[0])][self.edge_list.index(edge)] = -1 * \
                        float(component.parameter_list[0])
                elif component.edge_list.index(edge) == 1:
                    self.N_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = 1
                    self.N_matrix[self.edge_list.index(edge)][self.edge_list.index(component.edge_list[0])] = -1 * \
                        float(component.parameter_list[2])
                else:
                    pass
            elif component.type == "DBh":  # If component is Double Bipole with H prime Matrix
                if component.edge_list.index(edge) == 0:
                    self.N_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = 1
                    self.N_matrix[self.edge_list.index(edge)][self.edge_list.index(component.edge_list[1])] = -1 * \
                        float(component.parameter_list[1])
                elif component.edge_list.index(edge) == 1:
                    self.N_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = -1 * \
                        float(component.parameter_list[3])
                else:
                    pass
            elif component.type == "DBT":  # If component is Double Bipole with T Matrix
                if component.edge_list.index(edge) == 0:
                    self.N_matrix[self.edge_list.index(edge)][self.edge_list.index(edge.component.edge_list[1])] = \
                        float(component.parameter_list[2])
                elif component.edge_list.index(edge) == 1:
                    self.N_matrix[self.edge_list.index(edge)][self.edge_list.index(edge.component.edge_list[0])] = 1
                    self.N_matrix[self.edge_list.index(edge)][self.edge_list.index(edge)] = \
                        float(component.parameter_list[3])
                else:
                    pass
            else:
                pass

    def build_Z(self, h):
        for edge in self.edge_list:
            # WARNING: VERY COMPONENT-SPECIFIC
            component = edge.component
            if component.type == "R":  # If component is a Resistor
                pass
            elif component.type == "DCVS" or \
                 component.type == "DCCS":  # If component is a Direct Current Source
                self.Z_vector[self.edge_list.index(edge)] = float(component.parameter_list[0])
            elif component.type == "ACVS" or \
                 component.type == "ACCS":  # If component is an Alternating Current Source
                self.Z_vector[self.edge_list.index(edge)] = float(component.parameter_list[0]) *\
                                                            numpy.cos(float(component.parameter_list[2]))
            elif component.type == "SWVS" or \
                 component.type == "SWCS":  # If component is a Square Wave Source
                self.Z_vector[self.edge_list.index(edge)] = float(component.parameter_list[0])
            elif component.type == "TWVS" or \
                 component.type == "TWCS":  # If component is a Triangle Wave Source
                self.Z_vector[self.edge_list.index(edge)] = float(component.parameter_list[0])
            elif component.type == "L":  # If component is an Inductor
                if self.edge_list.index(edge) == self.edge_list.index(component.edge_list[0]):
                    pass
                else:
                    self.Z_vector[self.edge_list.index(edge)] = float(component.parameter_list[1])
            elif component.type == "C":  # If component is a Capacitor
                if self.edge_list.index(edge) == self.edge_list.index(component.edge_list[0]):
                    pass
                else:
                    self.Z_vector[self.edge_list.index(edge)] = float(component.parameter_list[1])
            elif component.type == "DBR":  # If component is Double Bipole with R Matrix
                if component.edge_list.index(edge) == 0:
                    self.Z_vector[self.edge_list.index(edge)] = float(component.parameter_list[4])
                elif component.edge_list.index(edge) == 1:
                    self.Z_vector[self.edge_list.index(edge)] = float(component.parameter_list[5])
                else:
                    pass
            elif component.type == "DBG":  # If component is Double Bipole with G Matrix
                if component.edge_list.index(edge) == 0:
                    self.Z_vector[self.edge_list.index(edge)] = float(component.parameter_list[4])
                elif component.edge_list.index(edge) == 1:
                    self.Z_vector[self.edge_list.index(edge)] = float(component.parameter_list[5])
                else:
                    pass
            elif component.type == "DBH":  # If component is Double Bipole with H Matrix
                if component.edge_list.index(edge) == 0:
                    self.Z_vector[self.edge_list.index(edge)] = float(component.parameter_list[4])
                elif component.edge_list.index(edge) == 1:
                    self.Z_vector[self.edge_list.index(edge)] = float(component.parameter_list[5])
                else:
                    pass
            elif component.type == "DBh":  # If component is Double Bipole with H prime Matrix
                if component.edge_list.index(edge) == 0:
                    self.Z_vector[self.edge_list.index(edge)] = float(component.parameter_list[4])
                elif component.edge_list.index(edge) == 1:
                    self.Z_vector[self.edge_list.index(edge)] = float(component.parameter_list[5])
                else:
                    pass
            elif component.type == "DBT":  # If component is Double Bipole with T Matrix
                if component.edge_list.index(edge) == 0:
                    self.Z_vector[self.edge_list.index(edge)] = float(component.parameter_list[4])
                elif component.edge_list.index(edge) == 1:
                    self.Z_vector[self.edge_list.index(edge)] = float(component.parameter_list[5])
                else:
                    pass
            else:
                pass


    def build_T(self):
        O_node_matrix = numpy.zeros((len(self.node_list) - 1, len(self.node_list) - 1))
        O_node_edge_matrix = numpy.zeros((len(self.node_list) - 1, len(self.edge_list)))
        O_edge_matrix = numpy.zeros((len(self.edge_list), len(self.edge_list)))
        I_edge_matrix = numpy.identity(len(self.edge_list))

        T_matrix_first = numpy.concatenate((O_node_matrix, O_node_edge_matrix), axis=1)
        T_matrix_first = numpy.concatenate((T_matrix_first, self.A_matrix), axis=1)
        T_matrix_second = numpy.concatenate((numpy.negative(self.A_matrix.T), I_edge_matrix), axis=1)
        T_matrix_second = numpy.concatenate((T_matrix_second, O_edge_matrix), axis=1)
        T_matrix_third = numpy.concatenate((O_node_edge_matrix.T, self.M_matrix), axis=1)
        T_matrix_third = numpy.concatenate((T_matrix_third, self.N_matrix), axis=1)

        T_matrix = numpy.concatenate((T_matrix_first, T_matrix_second), axis=0)
        return numpy.concatenate((T_matrix, T_matrix_third), axis=0)  # Return final T_matrix

    def solve(self):
        W_vector = numpy.concatenate((numpy.zeros(len(self.node_list) + len(self.edge_list) - 1), self.Z_vector))
        self.S_vector = numpy.linalg.solve(self.T_matrix, W_vector)
        return self.S_vector

    def print_edge_variables(self, solution_vec):
        offset = len(self.node_list) - 1
        for edge in self.edge_list:
            print("Voltage drop in " + edge.component.identifier + ", edge number " +\
                  str(edge.component.edge_list.index(edge)) + ": " + str(solution_vec[offset + self.edge_list.index(edge)]))
        for edge in self.edge_list:
            print("Current in " + edge.component.identifier + ", edge number " +
                  str(edge.component.edge_list.index(edge)) + ": " + str(solution_vec[offset + len(self.edge_list) + self.edge_list.index(edge)]))

    def update(self, solution_vec, time):  # Updates internal Z vector
        offset = len(self.node_list) - 1
        Z_vector_new = numpy.copy(self.Z_vector)
        for edge in self.edge_list:
            if edge.component.is_dynamic:  # Dynamic Components are tough, let's treat them separately
                if edge.component.type == "C":
                    if self.edge_list.index(edge) == self.edge_list.index(edge.component.edge_list[0]):
                        pass
                    else:
                        Z_vector_new[self.edge_list.index(edge.component.edge_list[1])] =\
                            solution_vec[offset + self.edge_list.index(edge)] + solution_vec[offset +
                            self.edge_list.index(edge.component.edge_list[0])]
                elif edge.component.type == "L":
                    if self.edge_list.index(edge) == self.edge_list.index(edge.component.edge_list[0]):
                        pass
                    else:
                        Z_vector_new[self.edge_list.index(edge.component.edge_list[1])] =\
                            solution_vec[offset + len(self.edge_list) + self.edge_list.index(edge)] + solution_vec[offset +
                            len(self.edge_list) + self.edge_list.index(edge.component.edge_list[0])]
                else:
                    pass
            else:  # These components may not be dynamic, but they still are time-variant
                if edge.component.type == "ACVS" or \
                   edge.component.type == "ACCS":  # Alternating Current
                    Z_vector_new[self.edge_list.index(edge)] = float(edge.component.parameter_list[0]) *\
                                                               numpy.cos(float(edge.component.parameter_list[1]) * time +
                                                                         float(edge.component.parameter_list[2]))
                if edge.component.type == "SWVS" or \
                   edge.component.type == "SWCS":  # Square Wave
                    if int(time / float(edge.component.parameter_list[2])) % 2 == 0:  # If we're in the first half period
                        Z_vector_new[self.edge_list.index(edge)] = float(edge.component.parameter_list[0])
                    else:  # Then we're in the second half period
                        Z_vector_new[self.edge_list.index(edge)] = float(edge.component.parameter_list[1])
                if edge.component.type == "TWVS" or \
                   edge.component.type == "TWCS":  # Triangle Wave
                    ratio = (float(edge.component.parameter_list[1]) - float(edge.component.parameter_list[0])) / \
                            float(edge.component.parameter_list[2])  # Slope of the triangle
                    delta = (time - float(edge.component.parameter_list[2]) * int(
                        time / float(edge.component.parameter_list[2])))  # Time passed since new half period
                    if int(time / float(edge.component.parameter_list[2])) % 2 == 0:  # If we're in the first half period
                        Z_vector_new[self.edge_list.index(edge)] = float(edge.component.parameter_list[0]) + delta * ratio
                    else:  # Then we're in the second half period
                        Z_vector_new[self.edge_list.index(edge)] = float(edge.component.parameter_list[1]) - delta * ratio
        self.Z_vector = Z_vector_new


class Simulation:
    def __init__(self, topology):
        self.topology = topology
        self.step = 0
        self.delta_time = 0
        self.time = 0
        self.set_time()
        self.topology.build_all(self.step)
        self.solution = []

    def set_time(self):
        print("Please choose a time interval (in seconds):")
        self.delta_time = float(input())
        print("Please choose an integration step:")
        self.step = float(input())

    def update_trackers(self):
        offset = len(self.topology.node_list) - 1
        for edge in self.topology.edge_list:
            if edge.component.is_tracked:
                if len(edge.component.edge_list) == 1:  # Bipoles are easy
                    edge.current_tracker.append(self.solution[offset + len(self.topology.edge_list) +
                                                                      self.topology.edge_list.index(edge)])
                    edge.voltage_tracker.append(self.solution[offset + self.topology.edge_list.index(edge)])
                else:
                    if edge.component.type == "C":  # If we're tracking a capacitor
                        if edge.component.edge_list.index(edge) == 0:  # Let's just consider the first edge
                            edge.current_tracker.append(self.solution[offset + len(self.topology.edge_list) +
                                                                      self.topology.edge_list.index(edge)])
                            edge.voltage_tracker.append(self.solution[offset + self.topology.edge_list.index(edge)] + \
                                                        self.solution[offset + self.topology.edge_list.index(
                                                            edge.component.edge_list[1])])
                    elif edge.component.type == "L":  # If we're tracking an inductor
                        if edge.component.edge_list.index(edge) == 0:  # Let's just consider the first edge
                            edge.current_tracker.append(self.solution[offset + len(self.topology.edge_list) +
                                                                      self.topology.edge_list.index(edge)] + \
                                                        self.solution[offset + len(self.topology.edge_list) +
                                                                      self.topology.edge_list.index(
                                                                          edge.component.edge_list[1])])
                            edge.voltage_tracker.append(self.solution[offset + self.topology.edge_list.index(edge)])
                    elif len(edge.component.edge_list) == 2 and len(edge.component.parameter_list) == 6:  # This should be all about double bipoles
                        # This happens once for every edge (two edges in double bipoles)
                        edge.current_tracker.append(self.solution[offset + len(self.topology.edge_list) +
                                                                      self.topology.edge_list.index(edge)])
                        edge.voltage_tracker.append(self.solution[offset + self.topology.edge_list.index(edge)])

    def plot_trackers(self):
        for edge in self.topology.edge_list:
            if edge.component.is_tracked:
                if len(edge.component.edge_list) == 2:
                    if edge.component.type == "C" or edge.component.type == "L":  # Now we're plotting a capacitor or an inductor
                        if edge.component.edge_list.index(edge) == 0:  # Consider the first edge so we don't do this twice
                            plot.figure(self.topology.edge_list.index(edge))

                            plot.subplot(211)
                            plot.plot(edge.current_tracker, 'y', label=edge.component.identifier + " current")
                            plot.plot(edge.voltage_tracker, 'b', label=edge.component.identifier + " voltage")
                            plot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                                       ncol=2, mode="expand", borderaxespad=0.)
                            plot.ylim([1.1 * min(min(edge.current_tracker), min(edge.voltage_tracker)),
                                       1.1 * max(max(edge.current_tracker), max(edge.voltage_tracker))])

                            plot.subplot(212)
                            plot.plot(edge.current_tracker, edge.voltage_tracker, 'r')

                            plot.show()
                    elif len(edge.component.parameter_list) == 6:  # This should be all about double bipoles
                        if edge.component.edge_list.index(edge) == 0:  # Consider the first edge so we don't do this twice
                            plot.figure(self.topology.edge_list.index(edge))

                            plot.subplot(211)
                            plot.plot(edge.current_tracker, 'y', label=edge.component.identifier + " current 1")
                            plot.plot(edge.voltage_tracker, 'b', label=edge.component.identifier + " voltage 1")
                            plot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                                       ncol=2, mode="expand", borderaxespad=0.)
                            plot.ylim([1.1 * min(min(edge.current_tracker), min(edge.voltage_tracker)),
                                       1.1 * max(max(edge.current_tracker), max(edge.voltage_tracker))])

                            plot.subplot(212)
                            plot.plot(edge.component.edge_list[1].current_tracker, 'y', label=edge.component.identifier + " current 2")
                            plot.plot(edge.component.edge_list[1].voltage_tracker, 'b', label=edge.component.identifier + " voltage 2")
                            plot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                                       ncol=2, mode="expand", borderaxespad=0.)
                            plot.ylim([1.1 * min(min(edge.component.edge_list[1].current_tracker), min(edge.component.edge_list[1].voltage_tracker)),
                                       1.1 * max(max(edge.component.edge_list[1].current_tracker), max(edge.component.edge_list[1].voltage_tracker))])

                            plot.show()
                    else:
                        pass
                elif len(edge.component.edge_list) == 1:  # Now we're plotting an easy 2-terminal, nothing to worry about
                    plot.figure(self.topology.edge_list.index(edge))
                    plot.plot(edge.current_tracker, 'y', label=edge.component.identifier + " current")
                    plot.plot(edge.voltage_tracker, 'b', label=edge.component.identifier + " voltage")
                    plot.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
                               ncol=2, mode="expand", borderaxespad=0.)
                    plot.ylim([1.1 * min(min(edge.current_tracker), min(edge.voltage_tracker)),
                               1.1 * max(max(edge.current_tracker), max(edge.voltage_tracker))])
                    plot.show()
                else:
                    pass

    def simulate(self):
        while self.time < self.delta_time:
            self.solution = self.topology.solve()  # At every instant, we solve the linear system
            self.time += self.step
            self.topology.update(self.solution, self.time)  # We then update the dynamic components
            self.update_trackers()

    def set_trackers(self, component_list):  # Indices of components in component_list
        print("Which components would you like to keep track of? (input netlist indices separated by blank spaces)")
        indices = input().split()
        for index in indices:
            component_list[int(index)].is_tracked = True
        if len(indices) != 0:
            return False
        else:
            return True


def checkNetList(netlist):  # check netlist status
    knownSymbols = ["DCVS", "DCCS", "SWVS", "SWCS", "ACVS", "ACCS", "TWVS", "TWCS"
                    "DBR", "DBG", "DBH", "DBh", "DBT",
                    "OP_AMP", "ID_T",
                    "VCVS", "VCCS", "CCVS", "CCCS",
                    "R", "L", "C",
                    "OC", "SC",
                    "GND"]
    for lines in netlist:
        line = lines.split()
        if line[0] not in knownSymbols:
            sys.exit("Error: unknown element found in netlist")

def init(filepath, COMPONENT_LIST, EDGE_LIST, NODE_LIST):
    with open(filepath, "r") as net:
        NETLIST_LINES = net.readlines()  # Read all lines and store them
    net.close()
    ground_node = "0"
    checkNetList(NETLIST_LINES)
    for line in NETLIST_LINES:  # Generate component list from netlist
        current_line = line.split()
        if len(current_line) > 2:
            current_nodes = current_line[2].split('#')
        current_params = []
        if len(current_line) > 3:
            current_params = current_line[3].split('#')
        # WE HAVE SOME NAMESAKES HERE
        if current_line[0] == "SC":  # Short Circuit Namesake
            COMPONENT_LIST.append(Component("DCVS", current_line[1], current_nodes, [], [0]))
        elif current_line[0] == "OC":  # Open Circuit Namesake
            COMPONENT_LIST.append(Component("DCCS", current_line[1], current_nodes, [], [0]))
        elif current_line[0] == "OP_AMP":  # Operational Amplifier Namesake
            COMPONENT_LIST.append(Component("DBT", current_line[1], current_nodes, [],
                                            [0, 0, 0, 0, 0, 0]))
        elif current_line[0] == "ID_T":  # Ideal Transformer Namesake
            COMPONENT_LIST.append(Component("DBT", current_line[1], current_nodes, [],
                                            [current_params[0], 0, 0, 1 / current_params[0], 0, 0]))
        elif current_line[0] == "VCVS":  # VCVS Namesake
            COMPONENT_LIST.append(Component("DBh", current_line[1], current_nodes, [],
                                            [0, 0, current_params[0], 0, 0, 0]))
        elif current_line[0] == "CCVS":  # CCVS Namesake
            COMPONENT_LIST.append(Component("DBR", current_line[1], current_nodes, [],
                                            [0, 0, current_params[0], 0, 0, 0]))
        elif current_line[0] == "VCCS":  # VCCS Namesake
            COMPONENT_LIST.append(Component("DBG", current_line[1], current_nodes, [],
                                            [0, 0, current_params[0], 0, 0, 0]))
        elif current_line[0] == "CCCS":  # CCCS Namesake
            COMPONENT_LIST.append(Component("DBH", current_line[1], current_nodes, [],
                                            [0, 0, current_params[0], 0, 0, 0]))
        elif current_line[0] == "GND":  # Ground Node Namesake
            ground_node = current_line[1]
        else:
            COMPONENT_LIST.append(Component(current_line[0], current_line[1], current_nodes, [], current_params))
        for node in current_nodes:
            NODE_LIST.append(node)

            # And now apply the ground node
    for node in NODE_LIST:
        if node == ground_node:
            NODE_LIST[NODE_LIST.index(node)] = "GND"
    for component in COMPONENT_LIST:
        for node in component.node_list:
            if node == ground_node:
                component.node_list[component.node_list.index(node)] = "GND"


            # COMPONENT-SPECIFIC EDGE GENERATION
    for component in COMPONENT_LIST:
        if component.type == "R" or \
                component.type == "DCVS" or \
                component.type == "DCCS" or \
                component.type == "ACVS" or \
                component.type == "ACCS" or \
                component.type == "SWVS" or \
                component.type == "SWCS" or \
                component.type == "TWVS" or \
                component.type == "TWCS":  # If component is "standard" bipole
            new_edge = Edge(component, component.node_list[0], component.node_list[1])
            component.edge_list.append(new_edge)
            EDGE_LIST.append(new_edge)
        elif component.type == "L":  # If component is an Inductor
            new_edge = Edge(component, component.node_list[0], component.node_list[1])
            component.edge_list.append(new_edge)
            EDGE_LIST.append(new_edge)
            new_edge = Edge(component, component.node_list[0], component.node_list[1])
            component.edge_list.append(new_edge)
            EDGE_LIST.append(new_edge)
            component.is_dynamic = True
        elif component.type == "C":  # If component is a Capacitor
            new_node = component.node_list[0] + "_th"
            component.node_list.append(new_node)
            NODE_LIST.append(new_node)
            new_edge = Edge(component, component.node_list[0], new_node)
            component.edge_list.append(new_edge)
            EDGE_LIST.append(new_edge)
            new_edge = Edge(component, new_node, component.node_list[1])
            component.edge_list.append(new_edge)
            EDGE_LIST.append(new_edge)
            component.is_dynamic = True
        elif component.type == "DBR" or \
             component.type == "DBG" or \
             component.type == "DBH" or \
             component.type == "DBh" or \
             component.type == "DBT":  # If component is Double Bipole with R Matrix
            new_edge = Edge(component, component.node_list[0], component.node_list[1])
            component.edge_list.append(new_edge)
            EDGE_LIST.append(new_edge)
            new_edge = Edge(component, component.node_list[2], component.node_list[3])
            component.edge_list.append(new_edge)
            EDGE_LIST.append(new_edge)
        else:
            pass

# ----------------------------------#
#           PROGRAM SETUP           #
# ----------------------------------#

COMPONENT_LIST = []
EDGE_LIST = []
NODE_LIST = []

print("Please insert Netlist filepath: ")
PATH = input()
init(PATH, COMPONENT_LIST, EDGE_LIST, NODE_LIST)

# ----------------------------------#
#           PROGRAM START           #
# ----------------------------------#

MAIN_CIRCUIT = Topology(NODE_LIST, EDGE_LIST)
MAIN_SIMULATION = Simulation(MAIN_CIRCUIT)

print("A matrix printout:")
print(MAIN_CIRCUIT.A_matrix)
print("\nM matrix printout:")
print(MAIN_CIRCUIT.M_matrix)
print("\nN matrix printout:")
print(MAIN_CIRCUIT.N_matrix)
print("\nT matrix printout:")
print(MAIN_CIRCUIT.T_matrix)
print("\nZ vector printout:")
print(MAIN_CIRCUIT.Z_vector)
print("\n\n\n\n")

if MAIN_SIMULATION.set_trackers(COMPONENT_LIST):
    MAIN_CIRCUIT.print_edge_variables(MAIN_CIRCUIT.solve())
else:
    MAIN_SIMULATION.simulate()
    MAIN_SIMULATION.plot_trackers()
