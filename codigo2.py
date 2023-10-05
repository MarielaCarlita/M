import copy
from itertools import combinations

vertices = [[1.0],[0.0],[0,1]]

simplice_0_dim = [0]
simplice_0_dim = [1,2]
simplice_0_dim = [0,1,2]

simplices = [[0],[1],[2],[0,1],[1,2]]

def dimen(simplex):
    return len(simplex)-1

def is_boundary(possible_boundary, simplex):
    if len(possible_boundary) +1 != len(simplex):
        return False
    extra = 0
    for k in range (len(possible_boundary)):
        if possible_boundary[k] == simplex[k+extra]:
            continue
        if extra == 0 and possible_boundary[k] == simplex[k+1]:
            extra = 1
            continue
        return False
    return True

def is_coboundary (possible_coboundary, simplex):
    return is_boundary(simplex, possible_coboundary)

filtracion_por_simplices = [
    [[1, 0], [0, 0], [0, 1]],
    [[0], [1], [2], [0, 1], [1, 2]]
]


#Para el acceso a los elementos de la ﬁltración creamos los siguientes métodos:

def get_vertices(filtration):
    return filtration[0]

def get_vertex(filtration, i):
    return filtration [0] [i]


#Para dibujar complejos simpliciales hasta en 3 dimensiones utilizaremos las siguientes funciones:

def plot_simplicial_complex_0_dim(simplicial_complex):
    number_of_vertices = len(get_vertices(simplicial_complex))
    simplex_0_dim = simplicial_complex[1][:number_of_vertices]
    return point([get_vertex(simplicial_complex, index[0]) for index in simplex_0_dim], color="red")

def plot_simplicial_complex_1_dim(simplicial_complex):
    number_of_vertices = len(get_vertices(simplicial_complex))
    simplex_1_dim = []
    for i in range(number_of_vertices, len(simplicial_complex[1])):
        if len(simplicial_complex[1][i]) != 2:
            break
        simplex_1_dim.append(simplicial_complex[1][i])
    return sum(line( [get_vertex(simplicial_complex, index[0]), get_vertex(simplicial_complex, index[1])], color="black"
) for index in simplex_1_dim)

def plot_simplicial_complex_2_dim(simplicial_complex):
    number_of_vertices = len(get_vertices(simplicial_complex))
    simplex_2_dim = []
    for i in range(number_of_vertices, len(simplicial_complex[1])):
        long = len(simplicial_complex[1][i])
        if long == 3:
            simplex_2_dim.append(simplicial_complex[1][i])
        elif long > 3:
            break
    return sum(polygon([get_vertex(simplicial_complex, index[0]),
                        get_vertex(simplicial_complex, index[1]),
                        get_vertex(simplicial_complex, index[2])])
            for index in simplex_2_dim)

def plot_simplicial_complex_3d (simplicial_complex):
    return sum([plot_simplicial_complex_0_dim(simplicial_complex),
                plot_simplicial_complex_1_dim(simplicial_complex),
                plot_simplicial_complex_2_dim(simplicial_complex)])


#Utilizamos la siguiente función para hacer print de matrices con una ﬁla por línea:

def matrix_print(cadena, M_t):
    print(cadena)
    for i in M_t:
        print(i)


def matrix_print_traspose(cadena, M_t):
    print(cadena)
    for i in traspose(M_t):
        print(i)

def get_boundary_matrix_t(rows, columns):
    return [[1 if is_boundary(rows[j], columns[i]) else 0
                    for j in range(len(rows))]
                    for i in range(len(columns))]

def get_coboundary_matrix_t(rows, columns):
    coboundary_t = []
    for i in range(len(columns)):
        coboundary_t.append([])
        for j in range(len(rows)):
            coboundary_t[i].append(1 if is_coboundary(rows[j], columns[i]) else 0)
    return coboundary_t

def fill_with_zeros(zeros_left, cadena, zeros_right):
    return [0 for _ in range(zeros_left)] + cadena + [0 for _ in range(zeros_right)]

def same_dimension(simplex_1, simplex_2):
    return len(simplex_1) == len(simplex_2)

def boundary_matrix_t(simplicial_complex):
    simplices_list = simplicial_complex[1]
    number_of_vertices = len(simplicial_complex[0])
    number_of_simplices = len(simplices_list)
    
    M_boundary_t = []
    empty_row = [0 for _ in range(number_of_simplices)]
    for _ in range(number_of_vertices):
        M_boundary_t.append(empty_row)

    start_row_index = 0
    start_column_index = number_of_vertices
    for column_index in range(start_column_index, number_of_simplices-1):
        if not same_dimension(simplices_list[column_index], simplices_list[column_index+1]):
            m_boundary_t = get_boundary_matrix_t(simplices_list[start_row_index: start_column_index], simplices_list[start_column_index:column_index+1])

#Fill with 0 boundary matrix
            for col_index in range(len(m_boundary_t)):
                M_boundary_t.append(fill_with_zeros(start_row_index,
                m_boundary_t[col_index], number_of_simplices - start_column_index))
            start_row_index = start_column_index
            start_column_index = column_index+1

    m_boundary_t = get_boundary_matrix_t(simplices_list[start_row_index: start_column_index], simplices_list[start_column_index:])

    #Fill with 0 boundary matrix

    for col_index in range(len(m_boundary_t)):
        M_boundary_t.append(fill_with_zeros(start_row_index, m_boundary_t[col_index], number_of_simplices - start_column_index))
        return M_boundary_t

def low(col):
    for possible_pivot in range(len(col)-1, -1, -1):
        if col[possible_pivot] == 1:
            return possible_pivot
    return -1

def lows(R_t):
    if len(R_t) == 0: return []
    pivots = []
    columns_to_search = [*range(len(R_t))]
    for row_index in range(len(R_t[0])-1, -1, -1):
        for i in range(len(columns_to_search)):
            if R_t[columns_to_search[i]][row_index] != 0:
                pivots.insert(0,row_index)
                columns_to_search.pop(i)
                break
    return pivots

def sum_cols_Z2(first_col, second_col):
    return [1 if first_col[i] != second_col[i] else 0
            for i in range(len(first_col))]

def create_V_matrix(dimen):
    return [[1 if i == j else 0 for j in range(dimen)] for i in range(dimen)]

def get_persistence_pairs_t(R_t):
    persistence_pairs = []
    for j in range(len(R_t)):
        lower = low(R_t[j])
        if lower != -1:
            persistence_pairs.append([lower, j])
    return persistence_pairs

def get_essential_indices_t(R_t, persistence_pairs):
    essential_indices = []
    birth = [pair[0] for pair in persistence_pairs]
    index = 0
    for essential in range(len(R_t)):
        if low(R_t[essential]) == -1 and essential not in birth:
            essential_indices.append(essential)
    return essential_indices

def divide_by_dimension(simplexwise_filtration):
    filtration_by_dimension = []
    p_simplices = []
    for i in range(len(simplexwise_filtration)-1):
        p_simplices.append(simplexwise_filtration[i].copy())
        if dimen(simplexwise_filtration[i]) != dimen(simplexwise_filtration[i+1]):
            filtration_by_dimension.append(p_simplices.copy())
            p_simplices = []
    p_simplices.append(simplexwise_filtration[-1].copy())
    filtration_by_dimension.append(p_simplices.copy())
    return filtration_by_dimension

def persistence_diagram_radious_points(number_of_simplices, diameters, persistence_pairs, essential_indices):
    points = []
    maximum = max(diameters)
    for pair in persistence_pairs:
        points.append([diameters[pair[0]], diameters[pair[1]]])

    for essential_index in essential_indices:
        points.append([diameters[essential_index], maximum])
    
    return points

def dimension_change_indices (split_filtration):
    indices = []
    summand = 0
    for d in split_filtration:
        indices.append(len(d)+summand)
        summand = indices[-1]
    return indices

def simplex_diameter (simplex, distance_m):
    if len(simplex) == 1:
        return 0
    return max([ distance_m[comb[0]][comb[-1]] for comb in combinations(simplex, 2) ])

def plot_persistence_diagram(simplexwise_filtration, persistence_pairs, essential_indices):
    split_filtration = divide_by_dimension(simplexwise_filtration)
    dimension_change_list = dimension_change_indices(split_filtration)
    
    number_of_pp = len(persistence_pairs)
    max_index = len(simplexwise_filtration)-1
    step = max_index/10
    
    points_by_dimension = [[] for _ in range(len(split_filtration))]
    for pair in persistence_pairs:
        for index in range(len(dimension_change_list)):
            if pair[0] < dimension_change_list[index]:
                points_by_dimension[index].append(pair)
                break
    for essential in essential_indices:
        for index in range(len(dimension_change_list)):
            if essential < dimension_change_list[index]:
                points_by_dimension[index].append(
                    [essential, max_index+(step/10)]
                )
                break

    colors = ["red", "blue", "green"]

    for i in range(3, len(points_by_dimension)):
        colors.append((random(),random(),random()))
    return sum(
        point(points_by_dimension[i], color=colors[i], size=15, legend_label=i)
        for i in range(len(points_by_dimension))
    )+line(
        [[0, max_index+(step/10)],[max_index, max_index+(step/10)]],
        color="black", linestyle='--',legend_label='infinity'
    )+line(
        [[0, 0],[max_index, max_index]],
        color="black", linestyle='--',legend_label='x=y'
    )

def count_and_remove_occurrences (element, element_list):
    count = 0
    index = 0
    while index < len(element_list):
        if element_list[index] == element:
            count += 1
            element_list.pop(index)
        else:
            index += 1
    return count

def test ():
    element = 2
    lista = [1,2,3,4,2]
    lista_exp = [1,3,4]
    count_exp = 2
    count = count_and_remove_occurrences (element, lista)
    assert lista == lista_exp
    assert count == count_exp

test()

def test ():
    element = 2
    lista = [1,2,2,4,2]
    lista_exp = [1,4]
    count_exp = 3
    count = count_and_remove_occurrences (element, lista)
    assert lista == lista_exp
    assert count == count_exp
test()

def plot_persistence_diagram_radious(simplexwise_filtration, distance_m, persistence_pairs, essential_indices):
    number_of_simplices = len(simplexwise_filtration)
    
    diameters = [simplex_diameter(simplex, distance_m) for simplex in simplexwise_filtration]
    
    max_diameter = max(diameters)
    step = max_diameter/10

    points = persistence_diagram_radious_points(
        number_of_simplices,
        diameters,
        persistence_pairs,
        essential_indices
    )

    number_of_pp = len(persistence_pairs)
    
    split_filtration = divide_by_dimension(simplexwise_filtration)
    dimension_change_list = dimension_change_indices(split_filtration)

    points_by_dimension = [[] for _ in range(len(split_filtration))]
    point_index = 0
    while point_index < number_of_pp:
        for index in range(len(dimension_change_list)):
            if persistence_pairs[point_index][0] < dimension_change_list[index]:
                points_by_dimension[index].append(points[point_index])
                break
        point_index += 1

    while point_index < len(points):
        for index in range(len(dimension_change_list)):
            if essential_indices[point_index - number_of_pp] < dimension_change_list[index]:
                points_by_dimension[index].append([points[point_index][0],points[point_index][1]+(step/10)])
                break
        point_index += 1

    colors = ["red", "blue", "green"]
    for i in range(3, len(points_by_dimension)):
        colors.append((random(),random(),random()))
    
    points_count = copy.deepcopy(points)
    multiplicity = []
    point_index = 0
    while point_index < len(points_count):
        p = points_count[point_index]
        count = count_and_remove_occurrences (p, points_count)
        if count != 1:
            multiplicity.append([count, [p[0],p[1]+(step/4)]])
    return line(
        [[0, max_diameter+(step/10)],[max_diameter, max_diameter+(step/10)]],
        color="black", linestyle='--',legend_label='infinity'
    )+line(
        [[0, 0],[max_diameter, max_diameter]],
        color="black", linestyle='--',legend_label='x=y'
    )+sum(
        point(points_by_dimension[i], color=colors[i], size=15, legend_label=i)
        for i in range(len(points_by_dimension))
    )+sum(
        text(multi[0], (multi[1]), fontsize='xx-small') for multi in multiplicity
    )

def traspose(matrix):
    if len(matrix) == 0: return []
    return [[matrix[i][j] for i in range(len(matrix))] for j in range(len(matrix[0]))]

def reduce_column_if_necessary(R_t, V_t, i):
    actual_low = low(R_t[i])
    j = 0
    while actual_low != -1 and j < i:
        if low(R_t[j]) == low(R_t[i]):
            R_t[i] = sum_cols_Z2 (R_t[i], R_t[j])
            V_t[i] = sum_cols_Z2 (V_t[i], V_t[j])
            actual_low = low(R_t[i])
            j = 0
        j += 1
    return R_t, V_t

def pHcol_Z2(M_boundary_t):
    assert len(M_boundary_t) > 0
    R_t = copy.deepcopy(M_boundary_t)
    number_of_columns = len(R_t)
    V_t = create_V_matrix(number_of_columns)
    for i in range(number_of_columns):
        R_t, V_t = reduce_column_if_necessary(R_t, V_t, i)
    return R_t, V_t

def indices_with_low_i (R_t, i):
    assert len(R_t) > 0
    indices = []
    for j in range(len(R_t)):
        if low(R_t[j]) == i:
            indices.append(j)
    return indices

def pHrow_Z2(M_boundary_t):
    assert len(M_boundary_t) > 0
    R_t = copy.deepcopy(M_boundary_t)
    dimen = len(R_t)
    V_t = create_V_matrix(dimen)
    for i in range(len(R_t[0])-1, -1, -1):
        indices = indices_with_low_i (R_t, i)
        if len(indices) > 0:
            reduced_col_index = indices[0]
            for j in indices[1:]:
                R_t[j] = sum_cols_Z2 (R_t[j], R_t[reduced_col_index])
                V_t[j] = sum_cols_Z2 (V_t[j], V_t[reduced_col_index])
    return R_t, V_t


cocycle = [0, [[1],[2]]]

cocycles = [[1,[[0,1]]], [1,[[1]]], [0,[[0],[1]]]]

def marked(cocycle):
    return cocycle[0] == 1

def unmarked(cocycle):
    return cocycle[0] == 0

def mark(cocycle):
    cocopy = copy.deepcopy(cocycle)
    cocopy[0] = 1
    return cocopy

def get_chain(cocycle):
    return cocycle[1]

#33 -> 61

def simplex_indices_with_this_coboundary (possible_coboundary, simplices_list):
    boundary_indices = []
    for j in range(len(simplices_list)):
        if is_coboundary(possible_coboundary, simplices_list[j]):
            boundary_indices.append(j)
    return boundary_indices

def is_coboundary_of_the_chain(possible_coboundary, simplices_cochain):
    indices = simplex_indices_with_this_coboundary(possible_coboundary, simplices_cochain)
    return len(indices) % 2 != 0

def cocycles_indices_with_this_coboundary(possible_coboundary, cocycles):
    cocycles_indices = []
    for j in range(len(cocycles)):
        if unmarked(cocycles[j]) and is_coboundary_of_the_chain(possible_coboundary, get_chain(cocycles[j])):
            cocycles_indices.append(j)
    return cocycles_indices

def pCoh_Z2(simplexwise_filtration):
    number_of_simplices = len(simplexwise_filtration)
    
    persistence_pairs=[]
    cocycles = []
    birth = []
    
    for i in range(number_of_simplices):
        new_simplex = simplexwise_filtration[i]

        former_cocycles_indices = cocycles_indices_with_this_coboundary(new_simplex, cocycles)

        if len(former_cocycles_indices) == 0:
            cocycles.insert(0,[0,[new_simplex]])
            birth.insert(0,i)
        else:
            p = former_cocycles_indices[0]
            for former_cocycle_index in former_cocycles_indices[1:]:
                cocycles[former_cocycle_index][1] = cocycles[former_cocycle_index][1] + cocycles[p][1]
            cocycles[p] = mark(cocycles[p])
            persistence_pairs.append([birth[p], i])
            cocycles.insert(0,[1,[new_simplex]])
            birth.insert(0,i)
    return persistence_pairs, cocycles, birth

def get_essential_indices_coh(cocycles, birth):
    essential_indices = []
    for i in range(len(cocycles)):
        if not marked(cocycles[i]):
            essential_indices.append(birth[i])
    return essential_indices

#35 -> 62 page
#pHcol con CLearing Column
def first_column_with_this_pivot(R_t, pivot_index):
    for j in range(len(R_t)):
        if low(R_t[j]) == pivot_index:
            return j
    return -1

def pHcol_Z2_cc(boundary_t):
    number_of_columns = len(boundary_t)
    R_t = copy.deepcopy(boundary_t)
    V_t = create_V_matrix(number_of_columns)
    persistence_pairs = []
    essential_indices = []

    for column_index in range(number_of_columns):
        pivot = low(R_t[column_index])
        reduced_column_index = first_column_with_this_pivot(R_t, pivot)
        while (
            pivot != -1 and
            reduced_column_index != column_index
        ):
            R_t[column_index] = sum_cols_Z2(R_t[column_index],R_t[reduced_column_index])

            V_t[column_index] = sum_cols_Z2(V_t[column_index],V_t[reduced_column_index])
            
            pivot = low(R_t[column_index])
            reduced_column_index = first_column_with_this_pivot(R_t, pivot)

        if pivot != -1:
            persistence_pairs.append([pivot, column_index])
        else:
            essential_indices.append(column_index)
    return V_t, R_t, persistence_pairs, essential_indices

def remove_elements_with_indices (array, indices_to_remove):
    new_array = []
    for i in range(len(array)):
        if i not in indices_to_remove:
            new_array.append(array[i])
    return new_array

def fix_pp_column_indices(persistence_pairs, pivots):
    fixed_pp = copy.deepcopy(persistence_pairs)
    for pair in fixed_pp:
        pivot_index = 0
        while pivot_index < len(pivots) and pair[1] >= pivots[pivot_index]:
            pivot_index += 1
            pair[1] += 1
    return fixed_pp

def fix_ess_column_indices(essential_indices, pivots):
    fixed_ess = copy.deepcopy(essential_indices)
    for ess_i in range(len(fixed_ess)):
        pivot_index = 0
        while pivot_index < len(pivots) and fixed_ess[ess_i] >= pivots[pivot_index]:
            pivot_index += 1
            fixed_ess[ess_i] += 1
    return fixed_ess

def add_empty_row_in_indices(matrix, indices):
    new_matrix = copy.deepcopy(matrix)
    for i in indices:
        new_matrix.insert(i,[0 for _ in new_matrix[0]])
    return new_matrix

def add_empty_col_in_indices(matrix, indices):
    new_matrix = traspose(matrix)
    for i in indices:
        new_matrix.insert(i,[0 for _ in new_matrix[0]])
    return traspose(new_matrix)

def extend_R_filling_in_zeros(R_t, pivots):
    pivots.sort()
    return add_empty_row_in_indices(R_t, pivots)

def extend_V_filling_in_zeros(V_t, pivots):
    new_V_t = add_empty_row_in_indices(V_t, pivots)
    new_V_t = add_empty_col_in_indices(new_V_t, pivots)
    return new_V_t

def clearing_columns (d_boundary_t, dp1_R_t, dp1_persistence_pairs):
    my_d_boundary_t = copy.deepcopy(d_boundary_t)

    pivots = lows(dp1_R_t)
    d_boundary_t_cleared = remove_elements_with_indices (my_d_boundary_t, pivots)

    d_V_t_cleared, d_R_t_cleared, d_persistence_pairs, d_essential_indices = pHcol_Z2_cc(d_boundary_t_cleared)
    d_persistence_pairs = fix_pp_column_indices(d_persistence_pairs, pivots)
    d_essential_indices = fix_ess_column_indices(d_essential_indices, pivots)
    d_R_t = extend_R_filling_in_zeros(d_R_t_cleared, pivots)

    if len(d_boundary_t[0]) != 0:
        d_V_t = extend_V_filling_in_zeros(d_V_t_cleared, pivots)

        for pair in dp1_persistence_pairs:
            d_V_t[pair[0]] = dp1_R_t[pair[1]]

        else:
            d_V_t = create_V_matrix(len(d_boundary_t))
    return d_V_t, d_R_t, d_persistence_pairs, d_essential_indices

def get_row_and_column_indices(filtration_by_dimension, columns_dimension):
    start_column_index = 0
    start_row_index = 0
    
    for d in range(columns_dimension-1):
        start_row_index += len(filtration_by_dimension[d])
    
    if columns_dimension != 0:
        start_column_index = start_row_index + len(filtration_by_dimension[columns_dimension-1])

    return start_row_index, start_column_index

def fix_pp (persistence_pairs, start_row_index, start_column_index):
    return [[pair[0]+start_row_index, pair[1]+start_column_index] for pair in persistence_pairs]

def fix_ess (essential_indices, start_column_index):
    return [essential_index+start_column_index for essential_index in essential_indices]

def pHcol_clearing_columns(simplicial_complex):
    simplexwise_filtration = copy.deepcopy(simplicial_complex[1])
    
    filtration_by_dimension = divide_by_dimension(simplexwise_filtration)
    
    total_simplices = len(simplexwise_filtration)
    number_of_dimensions = len(filtration_by_dimension)
    number_of_vertices = len(get_vertices(simplicial_complex))

    persistence_pairs = []
    essential_indices = []

    dp1_R_t = []
    dp1_pp = []

    for columns_dimension in range(number_of_dimensions-1, -1, -1):
        if columns_dimension == 0: rows = []

        else: rows = filtration_by_dimension[columns_dimension-1]
        columns = filtration_by_dimension[columns_dimension]
        d_boundary_t = get_boundary_matrix_t(rows, columns)
        d_V_t, d_R_t, d_pp, d_ess = clearing_columns(d_boundary_t, dp1_R_t, dp1_pp)

        start_row_index, start_column_index = get_row_and_column_indices(filtration_by_dimension, columns_dimension)
        persistence_pairs.extend(fix_pp(d_pp, start_row_index, start_column_index))

        essential_indices.extend(
            fix_ess (d_ess, start_column_index)
        )

        dp1_R_t = d_R_t
        dp1_pp = d_pp
    
    return persistence_pairs, essential_indices


#Intermedio Ripser

def first_column_with_this_pivot(R_t, pivot_index):
    for j in range(len(R_t)):
        if low(R_t[j]) == pivot_index:
            return j
    return -1

def pHcol_Z2_cc_cob(coboundary_t):
    number_of_columns = len(coboundary_t)

    R_t = copy.deepcopy(coboundary_t)
    V_t = create_V_matrix(number_of_columns)
    persistence_pairs = []
    birth_non_essential = []
    essential_indices = []

    for column_index in range(number_of_columns):
        pivot = low(R_t[column_index])
        while (
            pivot in birth_non_essential
        ):
            reduced_column_index = persistence_pairs[birth_non_essential.index(pivot)][1]
            R_t[column_index] = sum_cols_Z2(R_t[column_index],R_t[reduced_column_index])
            V_t[column_index] = sum_cols_Z2(V_t[column_index],V_t[reduced_column_index])

            pivot = low(R_t[column_index])

        if pivot != -1:
            persistence_pairs.append([pivot, column_index])
            birth_non_essential.append(pivot)
        else:
            essential_indices.append(column_index)

    return V_t, R_t, persistence_pairs, essential_indices

def clearing_columns_cob (d_boundary_t, dp1_R_t, dp1_persistence_pairs):
    my_d_boundary_t = copy.deepcopy(d_boundary_t)
    pivots = lows(dp1_R_t)
    d_boundary_t_cleared = remove_elements_with_indices (my_d_boundary_t, pivots)
    d_V_t_cleared, d_R_t_cleared, d_persistence_pairs, d_essential_indices = pHcol_Z2_cc_cob(d_boundary_t_cleared)

    d_persistence_pairs = fix_pp_column_indices(d_persistence_pairs, pivots)
    d_essential_indices = fix_ess_column_indices(d_essential_indices, pivots)
    d_R_t = extend_R_filling_in_zeros(d_R_t_cleared, pivots)

    if len(d_boundary_t[0]) != 0:
        d_V_t = extend_V_filling_in_zeros(d_V_t_cleared, pivots)
        for pair in dp1_persistence_pairs:
            d_V_t[pair[0]] = dp1_R_t[pair[1]]
    
    else:
        d_V_t = create_V_matrix(len(d_boundary_t))
    
    return d_V_t, d_R_t, d_persistence_pairs, d_essential_indices

def pHcol_clearing_columns_coboundary(simplicial_complex):
    simplexwise_filtration = copy.deepcopy(simplicial_complex[1])
    
    filtration_by_dimension = divide_by_dimension(simplexwise_filtration)
    
    total_simplices = len(simplexwise_filtration)
    number_of_dimensions = len(filtration_by_dimension)
    number_of_vertices = len(get_vertices(simplicial_complex))

    persistence_pairs = []
    essential_indices = []

    dp1_R_t = []
    dp1_pp = []

    for rows_dimension in range(1, number_of_dimensions):
        rows = copy.deepcopy(filtration_by_dimension[rows_dimension])
        number_of_rows = len(rows)
        
        columns = copy.deepcopy(filtration_by_dimension[rows_dimension-1])
        number_of_columns = len(columns)

        # Para la correcta computación de la coboundary
        rows.reverse()
        columns.reverse()
        d_coboundary_t = get_coboundary_matrix_t(rows, columns)
        d_V_t, d_R_t, d_pp, d_ess = clearing_columns_cob(d_coboundary_t, dp1_R_t, dp1_pp)

        d_persistence_pairs_cob = reverse_pairs_indices(number_of_rows -1, number_of_columns -1 , d_pp)

        d_essential_indices_cob = reverse_essential_indices(number_of_columns -1, d_ess)

        start_column_index, start_row_index = get_row_and_column_indices(filtration_by_dimension, rows_dimension)

        # Los pares de persistencia tienen que darse la vuelta
        d_persistence_pairs_cob = fix_pp(d_persistence_pairs_cob, start_row_index, start_column_index)
        d_persistence_pairs = reverse_each_pair(d_persistence_pairs_cob)

        persistence_pairs.extend(
            d_persistence_pairs
        )

        essential_indices.extend(
            fix_ess (d_essential_indices_cob, start_column_index)
        )

        dp1_R_t = d_R_t
        dp1_pp = d_pp
    
    # Rellenar los essential_indices que son todos los de mayor dimensión excepto los death
    birth_non_essential = [pair[1] for pair in d_persistence_pairs]
    for i in range(len(filtration_by_dimension[-1])):
        simplex_index = i + start_row_index
        if simplex_index not in birth_non_essential: essential_indices.append(simplex_index)
    return persistence_pairs, essential_indices

def reverse_pairs_indices (row_max_index, column_max_index, persistence_pairs):
    pp = []
    for pair in persistence_pairs:
        pair_reversed = [row_max_index - pair[0], column_max_index - pair[1]]
        pp.append(pair_reversed)
    return pp

def reverse_essential_indices (column_max_index, essential_indices):
    ess = []
    for essential_index in essential_indices:
        index_reversed = column_max_index - essential_index
        ess.append(index_reversed)
    return ess

def reverse_each_pair(persistence_pairs_cob):
    return [[pair[1],pair[0]] for pair in persistence_pairs_cob]

def get_pivot_and_compressed_column (column_simplex, rows, birth_non_essential):
    compressed_column = []
    check_if_is_pivot = True

    pivot = -1
    
    for row_index in range(len(rows)-1, -1, -1):
        if is_coboundary(rows[row_index], column_simplex):
            if check_if_is_pivot:
                check_if_is_pivot = False
                pivot = row_index
                if pivot not in birth_non_essential:
                    return pivot, []
            compressed_column.append(row_index)
        return pivot, compressed_column

def get_compressed_column (column_simplex, rows):
    compressed_column = []
    for row_index in range(len(rows)-1, -1, -1):
        if column_simplex[0] > rows[row_index][-1] or column_simplex[-1] < rows[row_index][0]:
            continue

        if is_coboundary(rows[row_index], column_simplex):
            compressed_column.append(row_index)
    return compressed_column

def get_pivot_from_compressed_column (compressed_column):
    if len(compressed_column) == 0:
        return -1

    possible_pivot = compressed_column[0]
    i = 1
    while i < len(compressed_column):
        if possible_pivot == -1:
            possible_pivot = compressed_column[i]

        elif compressed_column[i] != possible_pivot:
            return possible_pivot

        else:
            compressed_column.pop(0); compressed_column.pop(0)
            i -= 2
            possible_pivot = -1
        
        i += 1
    return possible_pivot

def sum_compressed_columns (compressed_column_1, compressed_column_2):
    result_compressed_column = []
    
    index_column_1 = 0
    index_column_2 = 0
    while index_column_1 < len(compressed_column_1) and index_column_2 < len(compressed_column_2):
        if compressed_column_1[index_column_1] >= compressed_column_2[index_column_2]:
            result_compressed_column.append(compressed_column_1[index_column_1])
            index_column_1 += 1

        else:
            result_compressed_column.append(compressed_column_2[index_column_2])
            index_column_2 += 1

    if index_column_1 < len(compressed_column_1):
        result_compressed_column.extend(compressed_column_1[index_column_1:])
    elif index_column_2 < len(compressed_column_2):
        result_compressed_column.extend(compressed_column_2[index_column_2:])

    return result_compressed_column

def compute_reduced_column(columns, rows, V_col):
    my_V_col = V_col.copy()
    column_index = get_pivot_from_compressed_column(my_V_col)
    if column_index == -1:
        return []
    
    column_index = get_pivot_from_compressed_column(my_V_col)
    reduced_column = get_compressed_column(columns[column_index], rows)

    number_of_pivots = 1
    column_index = get_pivot_from_compressed_column(my_V_col[number_of_pivots:])
    while column_index != -1:
        previously_sum_column = get_compressed_column(columns[column_index], rows)

        reduced_column = sum_compressed_columns(reduced_column, previously_sum_column)

        number_of_pivots +=1
        column_index = get_pivot_from_compressed_column(my_V_col[number_of_pivots:])
    return reduced_column

def sum_columns_to_reduce (rows, columns, column_to_reduce, V_reduced_column):
    working_column = column_to_reduce.copy()
    reduced_column = compute_reduced_column(columns, rows, V_reduced_column)
    working_column = sum_compressed_columns(working_column, reduced_column)
    return working_column

def create_compressed_V_matrix(number_of_columns):
    return [[i] for i in range(number_of_columns)]

def get_index_first_apparition (element, lista):
    index = 0
    while index < len(lista):
        if element == lista[index]:
            return index
        index += 1
    return -1

def pHcol_Z2_ripser(rows, columns):
    number_of_columns = len(columns)
    
    V_t = create_compressed_V_matrix(number_of_columns)
    persistence_pairs = []
    birth_non_essential = []
    essential_indices = []

    for column_index in range(number_of_columns):
        pivot, compressed_column = get_pivot_and_compressed_column (
            columns[column_index],
            rows,
            birth_non_essential
        )

        if pivot != -1 and len(compressed_column) == 0:
            persistence_pairs.append([pivot, column_index])
            birth_non_essential.append(pivot)

        else :
            birth_index = get_index_first_apparition(pivot, birth_non_essential)
            while (
                birth_index != -1 and pivot != -1
            ):
                reduced_column_index = persistence_pairs[birth_index][1]
                compressed_column = sum_columns_to_reduce (
                    rows, columns,
                    compressed_column,
                    V_t[reduced_column_index]
                )

            V_t[column_index] = sum_compressed_columns(
                V_t[column_index],
                V_t[reduced_column_index]
            )

            pivot = get_pivot_from_compressed_column (compressed_column)
            birth_index = get_index_first_apparition(pivot, birth_non_essential)

        if pivot != -1:
            persistence_pairs.append([pivot, column_index])
            birth_non_essential.append(pivot)
        else:
            essential_indices.append(column_index)
    return V_t, persistence_pairs, essential_indices

def get_pivots_sorted_from_pp (persistence_pairs):
    sorted_pivots = []
    for pair in persistence_pairs:
        bisect.insort(sorted_pivots, pair[0])
    return sorted_pivots

def clearing_columns_ripser (rows, columns, dp1_persistence_pairs):
    pivots = get_pivots_sorted_from_pp(dp1_persistence_pairs)
    columns_cleared = remove_elements_with_indices (columns, pivots)

    d_V_t_cleared, d_persistence_pairs, d_essential_indices = pHcol_Z2_ripser(rows, columns_cleared)
    d_persistence_pairs = fix_pp_column_indices(d_persistence_pairs, pivots)
    d_essential_indices = fix_ess_column_indices(d_essential_indices, pivots)
    return d_persistence_pairs, d_essential_indices

def ripser(simplicial_complex):
    simplexwise_filtration = copy.deepcopy(simplicial_complex[1])

    filtration_by_dimension = divide_by_dimension(simplexwise_filtration)

    total_simplices = len(simplexwise_filtration)
    number_of_dimensions = len(filtration_by_dimension)
    number_of_vertices = len(get_vertices(simplicial_complex))

    persistence_pairs = []
    essential_indices = []
    dp1_pp = []

    for rows_dimension in range(1, number_of_dimensions):
        rows = copy.deepcopy(filtration_by_dimension[rows_dimension])
        number_of_rows = len(rows)

        columns = copy.deepcopy(filtration_by_dimension[rows_dimension-1])
        number_of_columns = len(columns)

        # Para la correcta computación de la coboundary
        rows.reverse()
        columns.reverse()

        d_pp, d_ess = clearing_columns_ripser(rows, columns, dp1_pp)
        d_persistence_pairs_cob = reverse_pairs_indices(number_of_rows -1, number_of_columns -1 , d_pp)
        d_essential_indices_cob = reverse_essential_indices(number_of_columns-1, d_ess)

        start_column_index, start_row_index = get_row_and_column_indices(filtration_by_dimension, rows_dimension)

        # Los pares de persistencia tienen que darse la vuelta
        d_persistence_pairs_cob = fix_pp(d_persistence_pairs_cob, start_row_index, start_column_index)

        d_persistence_pairs = reverse_each_pair(d_persistence_pairs_cob)
        persistence_pairs.extend(
            d_persistence_pairs
        )

        essential_indices.extend(
            fix_ess (d_essential_indices_cob, start_column_index)
        )
        
        dp1_pp = d_pp

        # Rellenar los essential_indices que son todos los de mayor dimensión excepto los death
        birth_non_essential = [pair[1] for pair in d_persistence_pairs]
        for i in range(len(filtration_by_dimension[-1])):
            simplex_index = i + start_row_index
            if simplex_index not in birth_non_essential: essential_indices.append(simplex_index)
    
    return persistence_pairs, essential_indices


#Comparacion de Algoritmos

from scipy.spatial import distance_matrix
import bisect


def sort_by_diameter_and_index (chain_simplex, distance_m):
    def sorter (simplex):
        return (simplex_diameter (simplex, distance_m), simplex)
    
    return sorted(chain_simplex, key=sorter)

def is_valid_simplex(simplex, v, threshold, distance_m):
    for index in simplex:
        if distance_m[index][v] > threshold:
            return False
    return True

def add_next_dimension_simplices (last_sf, threshold, distance_m):
    d_simplexwise_filtration = []
    for v in range(len(distance_m)):
        for index in range(len(last_sf)):
            simplex = last_sf[index]
            if (
                v not in simplex and
                is_valid_simplex(simplex, v, threshold, distance_m)
            ):
                new_simplex = last_sf[index].copy()
                if (
                    v not in simplex and
                    is_valid_simplex(simplex, v, threshold, distance_m)
                ):
                    new_simplex = last_sf[index].copy()
                    bisect.insort(new_simplex, v)

                    if new_simplex not in d_simplexwise_filtration:
                        d_simplexwise_filtration.append(new_simplex)
    return d_simplexwise_filtration

def create_rips_simplexwise_filtration(vertices, threshold, max_dimension, distance_m):
    simplexwise_filtration = [[i] for i in range(len(vertices))]
    
    last_simplexwise_filtration = copy.deepcopy(simplexwise_filtration)
    for _ in range(1, max_dimension):
        d_simplexwise_filtration =add_next_dimension_simplices(last_simplexwise_filtration, threshold, distance_m)
        
        if len(d_simplexwise_filtration) == 0:
            break

        d_simplexwise_filtration = sort_by_diameter_and_index(d_simplexwise_filtration, distance_m)
        last_simplexwise_filtration = copy.deepcopy(d_simplexwise_filtration)
        simplexwise_filtration.extend(copy.deepcopy(d_simplexwise_filtration))
    return simplexwise_filtration

from time import time

#69 -> page 76
def print_times(matrix_times, complexes_names, alg_names):
    number_of_complexes = len(complexes_names)
    number_of_algorithms = len(alg_names)

    assert len(matrix_times) == number_of_algorithms
    assert len(matrix_times[0]) == number_of_complexes

    to_plot = []
    for alg_index in range(number_of_algorithms):
        alg_line = []
        for sc_index in range(number_of_complexes):
            alg_line.append([sc_index+1, matrix_times[alg_index][sc_index]])

        new_color = (random(), random(), random())
        to_plot.append(line(alg_line, color=new_color))
        to_plot.append(
            text(
                alg_names[alg_index],
                ([number_of_complexes+1, matrix_times[alg_index][-1]]),
                color=new_color
            )
        )
    maxi = 0
    for row in matrix_times:
        aux = max(row)
        maxi = aux if aux > maxi else maxi

    for sc_index in range(number_of_complexes):
        to_plot.append(
            text(
                complexes_names[sc_index],
                ([sc_index+1, maxi])
            )
        )

    return sum(to_plot)

def is_simplex_in_pp(simplex_index, persistence_pairs):
    for pair in persistence_pairs:
        if simplex_index in pair:
            return True
    return False

def each_simplex_in_pp_or_ess(total_simplices, persistence_pairs, essential_indices):
    for simplex in range(total_simplices):
        if simplex not in essential_indices and not is_simplex_in_pp(simplex,persistence_pairs):
            return False
    return True

def comprobar_persistencia (persistence_pairs_returned, essential_indices_returned):
    for pp_index in range(1, len(persistence_pairs_returned)):
        assert len(persistence_pairs_returned[pp_index-1]) == len(persistence_pairs_returned[pp_index])

    for ess_index in range(1, len(essential_indices_returned)):
        assert len(essential_indices_returned[ess_index-1]) == len(essential_indices_returned[ess_index])

    for ess_index in range(1, len(persistence_pairs_returned)):
        for pair in persistence_pairs_returned[0]:
            assert pair in persistence_pairs_returned[pp_index]

    for ess_index in range(1, len(essential_indices_returned)):
        for index in essential_indices_returned[0]:
            assert essential in essential_indices_returned[ess_index]

    total_simplices = 2*len(persistence_pairs_returned[0]) + len(essential_indices_returned[0])
    assert each_simplex_in_pp_or_ess(
        total_simplices,
        persistence_pairs_returned[0],
        essential_indices_returned[0]
    )
 
 # devolvera una lista con los pares de persistencia
PHCOL = 1
PHROW = 2
PHCOL_CC = 3
PHCOL_CC_COH = 4
PCOH =5
RIPSER = 6

def measure_times(simplicial_complex, algorithms_list):

    persistence_pairs_all = []

    essential_indices_all = []

    times = []

    if PHCOL in algorithms_list:
        t0 = time()

        M_boundary_t = boundary_matrix_t(simplicial_complex)
        R_t, V_t = pHcol_Z2(M_boundary_t)
        persistence_pairs = get_persistence_pairs_t(R_t)
        essential_indices = get_essential_indices_t(R_t, persistence_pairs)

        t1 = time()

        times.append(t1-t0)
        persistence_pairs_all.append(persistence_pairs)
        essential_indices_all.append(essential_indices)

    if PHROW in algorithms_list:
        t0 = time()

        M_boundary_t = boundary_matrix_t(simplicial_complex)
        R_t, V_t = pHrow_Z2(M_boundary_t)
        persistence_pairs = get_persistence_pairs_t(R_t)
        essential_indices = get_essential_indices_t(R_t, persistence_pairs)

        t1 = time()

        times.append(t1-t0)
        persistence_pairs_all.append(persistence_pairs)
        essential_indices_all.append(essential_indices)

    if PHCOL_CC in algorithms_list:
        t0 = time()

        persistence_pairs, essential_indices = pHcol_clearing_columns(simplicial_complex)

        t1 = time()

        times.append(t1-t0)
        persistence_pairs_all.append(persistence_pairs) 
        essential_indices_all.append(essential_indices)

    if PHCOL_CC_COH in algorithms_list:
        t0 = time()

        persistence_pairs, essential_indices = pHcol_clearing_columns_coboundary(simplicial_complex)

        t1 = time()

        times.append(t1-t0)
        persistence_pairs_all.append(persistence_pairs) 
        essential_indices_all.append(essential_indices)

    if PCOH in algorithms_list:
        t0 = time()

        simplexwise_filtration = simplicial_complex[1]
        persistence_pairs, cocycles, birth = pCoh_Z2(simplexwise_filtration)
        essential_indices = get_essential_indices_coh(cocycles, birth)

        t1 = time()

        times.append(t1-t0)
        persistence_pairs_all.append(persistence_pairs)
        essential_indices_all.append(essential_indices)

    if RIPSER in algorithms_list:
        t0 = time()

        persistence_pairs, essential_indices = ripser(simplicial_complex)

        times.append(t1-t0)

        persistence_pairs_all.append(persistence_pairs)
        essential_indices_all.append(essential_indices)

    return persistence_pairs_all, essential_indices_all, times

# Matriz que almacenara los tiempos
matrix_times = []
algorithms = [PHCOL, PHROW, PHCOL_CC, PHCOL_CC_COH, PCOH, RIPSER]
T1 = [0, 0, 0]; T2 = [2, 0 ,0]; T3 = [1, 0, 1.732]; T4 = [0.866, 1.732,0,866];

vertices = [T1, T2, T3, T4]

distance_m = distance_matrix(vertices, vertices)
simplexwise_filtration = create_rips_simplexwise_filtration(vertices, 3, 10, distance_m)

simplicial_complex_1 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_1)
plot_simplicial_complex_3d(simplicial_complex)

persistence_pairs_all, essential_indices_all, times_sc = measure_times(simplicial_complex, algorithms)
comprobar_persistencia(persistence_pairs_all, essential_indices_all)
print(persistence_pairs_all[0])
# cometario?? [[3,4], [2,5], [1,6], [1,10], [8,11], [9,12], [13,14]]
plot_persistence_diagram_radious(
    simplexwise_filtration,
    distance_m,
    persistence_pairs_all[0],
    essential_indices_all[0]
)
matrix_times.append(times_sc.copy())

#segundo ejemplo de complejo simplicial
CS0 = [0,2]; CS1 = [1.45,3.7]; CS2 = [1.4, -0.5]; CS3 = [2.3,1.6]; CS4 = [3.6,-0.1];
CS5 = [3.9,3]; CS6 = [5.1,1]; CS7 = [6.4, -0.5]; CS8 = [6.5,4.3]; CS9 = [7.2, 1.6];

threshold = 2.7
vertices = [CS0, CS1, CS2, CS3, CS4, CS5, CS6, CS7, CS8, CS9]

distance_m = distance_matrix(vertices, vertices)
simplexwise_filtration = create_rips_simplexwise_filtration(vertices, threshold, 10, distance_m)
simplicial_complex_2 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_2)
plot_simplicial_complex_3d(simplicial_complex)

persistence_pairs_all, essential_indices_all, times_sc = measure_times(
    simplicial_complex, 
    algorithms
)
comprobar_persistencia(persistence_pairs_all, essential_indices_all)
plot_persistence_diagram_radious(
    simplexwise_filtration,
    distance_m,
    persistence_pairs_all[0],
    essential_indices_all[0]
)
matrix_times.append(times_sc.copy())

CS0 = [-2, 0, 0]; CS1 = [0, 0, 0]; CS2 = [2, 0, 0]; CS3 = [3, 2, 0]; CS4 = [4, 0, 0];
CS5 = [6, 0.3, -1]; CS6 = [6.1, 1.9 , 0]; CS7 = [7, 0, 0]; CS8 = [7.8, 0.6, -1]; 

vertices = [CS0, CS1, CS2, CS3, CS4, CS5, CS6, CS7, CS8]
chain_simplex_dimen_0 = [[0],[1],[2],[3],[4],[5],[6],[7],[8]]
chain_simplex_dimen_1 = [[0,1],[1,2],[1,3],[2,3],[2,4],[3,4],[3,6],[4,7],[5,6],[5,7],[5,8],[6,7],[6,8],[7,8]]
chain_simplex_dimen_2 = [[2,3,4],[5,6,7],[5,6,8],[5,7,8],[6,7,8]]
chain_simplex_dimen_3 = [[5,6,7,8]]
simplexwise_filtration = chain_simplex_dimen_0 + chain_simplex_dimen_1 + chain_simplex_dimen_2 + chain_simplex_dimen_3
simplicial_complex_3 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_3)
plot_simplicial_complex_3d(simplicial_complex)

persistence_pairs_all, essential_indices_all, times_sc = measure_times(
    simplicial_complex, 
    algorithms
)
comprobar_persistencia(persistence_pairs_all, essential_indices_all)
plot_persistence_diagram(simplexwise_filtration, persistence_pairs_all[0], essential_indices_all[0])
matrix_times.append(times_sc.copy())
print_times(traspose(matrix_times), ["SC_1", "SC_2", "SC_3"], ["PHCOL", "PHROW", "PHCOL_CC", "PHCOL_CC_COH", "PCOH", "RIPSER"])

#CREANDO COMPLEJOS SIMPLICIALES CON UN NUMERO PEQUEÑO DE COMPONENTES CONEXAS
matrix_times = []
algorithms = [PHCOL, PHROW, PHCOL_CC, PHCOL_CC_COH, PCOH, RIPSER]
vertices_step = 1
number_of_vertices = 10
threshold = 10
max_simplex_dimension = 10

vertices = [[random(), random()] for _ in range(number_of_vertices)]

distance_m = distance_matrix(vertices, vertices)
simplexwise_filtration = create_rips_simplexwise_filtration(vertices, threshold, max_simplex_dimension, distance_m)

simplicial_complex_1 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_1)
plot_simplicial_complex_3d(simplicial_complex)  

persistence_pairs_all, essential_indices_all, times_sc = measure_times(
    simplicial_complex, 
    algorithms
)
comprobar_persistencia(persistence_pairs_all, essential_indices_all)
plot_persistence_diagram_radious(
    simplexwise_filtration,
    distance_m,
    persistence_pairs_all[0],
    essential_indices_all[0]
)

matrix_times.append(times_sc)
number_of_vertices = number_of_vertices + vertices_step
threshold = 10
max_simplex_dimension = 3

vertices = [[random(), random()] for _ in range(number_of_vertices)]

distance_m = distance_matrix(vertices, vertices)
simplexwise_filtration = create_rips_simplexwise_filtration(vertices, threshold, max_simplex_dimension, distance_m)

simplicial_complex_2 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_2)
plot_simplicial_complex_3d(simplicial_complex)

persistence_pairs_all, essential_indices_all, times_sc = measure_times(
    simplicial_complex, 
    algorithms
)
comprobar_persistencia(persistence_pairs_all, essential_indices_all)
plot_persistence_diagram_radious(
    simplexwise_filtration,
    distance_m,
    persistence_pairs_all[0],
    essential_indices_all[0]
)

matrix_times.append(times_sc)
number_of_vertices = number_of_vertices + vertices_step
threshold = 10
max_simplex_dimension = 3

vertices = [[random(), random()] for _ in range(number_of_vertices)]

distance_m = distance_matrix(vertices, vertices)
simplexwise_filtration = create_rips_simplexwise_filtration(vertices, threshold, max_simplex_dimension, distance_m)

simplicial_complex_3 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_3)
plot_simplicial_complex_3d(simplicial_complex)

persistence_pairs_all, essential_indices_all, times_sc = measure_times(
    simplicial_complex, 
    algorithms
)
comprobar_persistencia(persistence_pairs_all, essential_indices_all)
plot_persistence_diagram_radious(
    simplexwise_filtration,
    distance_m,
    persistence_pairs_all[0],
    essential_indices_all[0]
)

matrix_times.append(times_sc)
number_of_vertices = number_of_vertices + vertices_step
threshold = 10
max_simplex_dimension = 3

vertices = [[random(), random()] for _ in range(number_of_vertices)]

distance_m = distance_matrix(vertices, vertices)
simplexwise_filtration = create_rips_simplexwise_filtration(vertices, threshold, max_simplex_dimension, distance_m)

simplicial_complex_1 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_1)
plot_simplicial_complex_3d(simplicial_complex)

persistence_pairs_all, essential_indices_all, times_sc = measure_times(
    simplicial_complex, 
    algorithms
)
comprobar_persistencia(persistence_pairs_all, essential_indices_all)
plot_persistence_diagram_radious(
    simplexwise_filtration,
    distance_m,
    persistence_pairs_all[0],
    essential_indices_all[0]
)

matrix_times.append(times_sc)
number_of_vertices = number_of_vertices + vertices_step
threshold = 10
max_simplex_dimension = 3

vertices = [[random(), random()] for _ in range(number_of_vertices)]

distance_m = distance_matrix(vertices, vertices)
simplexwise_filtration = create_rips_simplexwise_filtration(vertices, threshold, max_simplex_dimension, distance_m)

simplicial_complex_1 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_1)
plot_simplicial_complex_3d(simplicial_complex)

persistence_pairs_all, essential_indices_all, times_sc = measure_times(
    simplicial_complex, 
    algorithms
)
comprobar_persistencia(persistence_pairs_all, essential_indices_all)
plot_persistence_diagram_radious(
    simplexwise_filtration,
    distance_m,
    persistence_pairs_all[0],
    essential_indices_all[0]
)

matrix_times.append(times_sc)
print_times(traspose(matrix_times)
            ["sk_10_3", "sk_11_3", "sk_12_3", "sk_14_3"], 
            ["PHCOL", "PHROW", "PHCOL_CC", "PHCOL_CC_COH", "PCOH", "RIPSER"])

matrix_times = []
algorithms = [PHCOL_CC, PHCOL_CC_COH, PCOH, RIPSER]
threshold = 10
max_simplex_dimension = 3
vertices_step = 2
number_of_vertices = 12

vertices = [[random(), random()] for _ in range(number_of_vertices)]

distance_m = distance_matrix(vertices, vertices)
simplexwise_filtration = create_rips_simplexwise_filtration(vertices, threshold, max_simplex_dimension, distance_m)

simplicial_complex_1 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_1)
plot_simplicial_complex_3d(simplicial_complex)

persistence_pairs_all, essential_indices_all, times_sc = measure_times(
    simplicial_complex, 
    algorithms
)
comprobar_persistencia(persistence_pairs_all, essential_indices_all)
plot_persistence_diagram_radious(
    simplexwise_filtration,
    distance_m,
    persistence_pairs_all[0],
    essential_indices_all[0]
)

matrix_times.append(times_sc)
number_of_vertices = number_of_vertices + vertices_step

vertices = [[random(), random()] for _ in range(number_of_vertices)]

distance_m = distance_matrix(vertices, vertices)
simplexwise_filtration = create_rips_simplexwise_filtration(vertices, threshold, max_simplex_dimension, distance_m)

simplicial_complex_2 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_2)
plot_simplicial_complex_3d(simplicial_complex)

persistence_pairs_all, essential_indices_all, times_sc = measure_times(
    simplicial_complex, 
    algorithms
)
comprobar_persistencia(persistence_pairs_all, essential_indices_all)
plot_persistence_diagram_radious(
    simplexwise_filtration,
    distance_m,
    persistence_pairs_all[0],
    essential_indices_all[0]
)

#codigo hasta  la pag 100
matrix_times.append(times_sc)
number_of_vertices = number_of_vertices + vertices_step 

vertices = [[random(), random()] for _ in range(number_of_vertices)]

distance_m = distance_matrix(vertices, vertices)
simplexwise_filtration = create_rips_simplexwise_filtration(vertices, max_simplex_dimension, distance_m)

simplicial_complex_3 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_3)
plot_simplicial_complex_3d(simplicial_complex)

persistence_pairs_all, essential_indices_all, times_sc = measure_times(
    simplicial_complex,
    algorithms
)

comprobar_persistencia(persistence_pairs_all, essential_indices_all)
plot_persistence_diagram_radious(
    simplexwise_filtration,
    distance_m,
    persistence_pairs_all[0],
    essential_indices_all[0]
)

matrix_times.append(times_sc)
number_of_vertices = number_of_vertices + vertices_step

vertices = [[random(), random()] for _ in range(number_of_vertices)]

distance_m = distance_matrix(vertices, vertices)
simplexwise_filtration = create_rips_simplexwise_filtration(vertices, threshold,
    max_simplex_dimension, distance_m)

simplicial_complex_1 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_1)
plot_simplicial_complex_3d(simplicial_complex)

persistence_pairs_all, essential_indices_all, times_sc = measure_times(
    simplicial_complex,
    algorithms
)

comprobar_persistencia(persistence_pairs_all, essential_indices_all)
plot_persistence_diagram_radious(
    simplexwise_filtration,
    distance_m,
    persistence_pairs_all[0],
    essential_indices_all[0]
)

matrix_times.append(times_sc)
number_of_vertices = number_of_vertices + vertices_step

vertices = [[random(), random()] for _ in range(number_of_vertices)]

distance_m = distance_matrix(vertices, vertices)
simplexwise_filtration = create_rips_simplexwise_filtration(
    vertices, threshold, max_simplex_dimension, distance_m)

simplicial_complex_1 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_1)
plot_simplicial_complex_3d(simplicial_complex)

persistence_pairs_all, essential_indices_all, times_sc = measure_times(
    simplicial_complex,
    algorithms
)

comprobar_persistencia(persistence_pairs_all, essential_indices_all)
plot_persistence_diagram_radious(
    simplexwise_filtration,
    distance_m,
    persistence_pairs_all[0],
    essential_indices_all[0]
)

matrix_times.append(times_sc)
print_times(traspose(matrix_times),
            ["sk_12_3", "sk_14_3", "sk_16_3", "sk_18_3", "sk_20_3"],
            ["PHCOL_CC", "PHCOL_CC_COH", "PCOH", "RIPSER"])

#dejamos de utilizar el algoridmo PHCOL_CC
matrix_times = []
algorithms = [ PHCOL_CC_COH, PCOH, RIPSER]

threshold = 10
max_simplex_dimension = 3
vertices_step = 5

number_of_vertices = 15

vertices = [[random(), random()] for _ in range(number_of_vertices)]

distance_m = distance_matrix(vertices, vertices)
simplexwise_filtration = create_rips_simplexwise_filtration(
    vertices, threshold, max_simplex_dimension, distance_m)

simplicial_complex_1 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_1)

plot_simplicial_complex_3d(simplicial_complex)

persistence_pairs_all, essential_indices_all, times_sc = measure_times(
    simplicial_complex,
    algorithms
)

comprobar_persistencia(persistence_pairs_all, essential_indices_all)
plot_persistence_diagram_radious(
    simplexwise_filtration, 
    distance_m,
    persistence_pairs_all[0], 
    essential_indices_all[0]
)

matrix_times.append(times_sc)
number_of_vertices = number_of_vertices + vertices_step

vertices = [[random(), random()] for _ in range(number_of_vertices)]

distance_m = distance_matrix(vertices, vertices)
simplexwise_filtration = create_rips_simplexwise_filtration(
    vertices, threshold, max_simplex_dimension, distance_m)

simplicial_complex_2 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_2)
plot_simplicial_complex_3d(simplicial_complex)

persistence_pairs_all, essential_indices_all, times_sc = measure_times(
    simplicial_complex,
    algorithms
)

comprobar_persistencia(persistence_pairs_all, essential_indices_all)
plot_persistence_diagram_radious(
    simplexwise_filtration, 
    distance_m,
    persistence_pairs_all[0], 
    essential_indices_all[0]
)

matrix_times.append(times_sc)
number_of_vertices = number_of_vertices + vertices_step

vertices = [[random(), random()] for _ in range(number_of_vertices)]

distance_m = distance_matrix(vertices, vertices)
simplexwise_filtration = create_rips_simplexwise_filtration(
    vertices, threshold, max_simplex_dimension, distance_m)

simplicial_complex_3 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_3)
plot_simplicial_complex_3d(simplicial_complex)

persistence_pairs_all, essential_indices_all, times_sc = measure_times(
    simplicial_complex,
    algorithms
)

comprobar_persistencia(persistence_pairs_all, essential_indices_all)
plot_persistence_diagram_radious(
    simplexwise_filtration, 
    distance_m,
    persistence_pairs_all[0], 
    essential_indices_all[0]
)

matrix_times.append(times_sc)
number_of_vertices = number_of_vertices + vertices_step

vertices = [[random(), random()] for _ in range(number_of_vertices)]

distance_m = distance_matrix(vertices, vertices)
simplexwise_filtration = create_rips_simplexwise_filtration(
    vertices, threshold, max_simplex_dimension, distance_m)

simplicial_complex_1 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_1)
plot_simplicial_complex_3d(simplicial_complex)

persistence_pairs_all, essential_indices_all, times_sc = measure_times(
    simplicial_complex,
    algorithms
)

comprobar_persistencia(persistence_pairs_all, essential_indices_all)
plot_persistence_diagram_radious(
    simplexwise_filtration, 
    distance_m,
    persistence_pairs_all[0], 
    essential_indices_all[0]
)

matrix_times.append(times_sc)
number_of_vertices = number_of_vertices + vertices_step

vertices = [[random(), random()] for _ in range(number_of_vertices)]

distance_m = distance_matrix(vertices, vertices)
simplexwise_filtration = create_rips_simplexwise_filtration(
    vertices, threshold, max_simplex_dimension, distance_m)

simplicial_complex_1 = [copy.deepcopy(vertices), copy.deepcopy(simplexwise_filtration)]
simplicial_complex = copy.deepcopy(simplicial_complex_1)
plot_simplicial_complex_3d(simplicial_complex)

persistence_pairs_all, essential_indices_all, times_sc = measure_times(
    simplicial_complex,
    algorithms
)

comprobar_persistencia(persistence_pairs_all, essential_indices_all)
plot_persistence_diagram_radious(
    simplexwise_filtration, 
    distance_m,
    persistence_pairs_all[0], 
    essential_indices_all[0]
)

matrix_times.append(times_sc)
print_times(traspose(matrix_times),
            ["sk_15_3", "sk_20_3", "sk_25_3", "sk_30_3", "sk_35_3"],
            ["PHCOL_CC_COH", "PCOH", "RIPSER"])
# fin
