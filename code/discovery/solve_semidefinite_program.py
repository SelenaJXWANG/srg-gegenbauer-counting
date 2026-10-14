"""Fresh cone SDP solve of the exact article model in data.json.

No historical certificate or weight direction is loaded.  This program
records numerical exploration only; exact reconstruction is a separate step.
"""
from pathlib import Path
from fractions import Fraction
import argparse
import hashlib
import json
import sys
import time

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'python'))
import cvxpy as cp
import numpy as np


def floating(value):
    if isinstance(value, list):
        return [floating(x) for x in value]
    return float(Fraction(value))


def statistics(problem):
    s = problem.solver_stats
    return {'solver_name': s.solver_name, 'solve_time': s.solve_time,
            'setup_time': s.setup_time, 'num_iters': s.num_iters,
            'extra_stats': repr(s.extra_stats)}


def run(data_path, output_path, solver='CLARABEL', variable_scale=10**10):
    started = time.time()
    data_bytes = data_path.read_bytes()
    data = json.loads(data_bytes)
    affine = np.array(floating(data['affine']))
    matrices = [np.array(floating(b['M'])) for b in data['blocks']]
    coordinate_scale = np.array([1.] + [float(variable_scale)] * (affine.shape[1]-1))
    scaled_affine = affine * coordinate_scale
    row_scales = np.max(np.abs(scaled_affine), axis=1)
    assert np.all(row_scales > 0)
    scaled_affine = scaled_affine / row_scales[:, None]
    block_scales = [float(np.max(np.abs(M * coordinate_scale[:, None, None])))
                    for M in matrices]
    assert min(block_scales) > 0
    scaled_matrices = [M * coordinate_scale[:, None, None] / scale
                       for M, scale in zip(matrices, block_scales)]
    q = affine.shape[1] - 1
    options = ({'tol_gap_abs': 1e-10, 'tol_feas': 1e-10,
                'tol_gap_rel': 1e-10, 'max_iter': 500}
               if solver == 'CLARABEL' else {'eps': 1e-7, 'max_iters': 250000})
    t = cp.Variable(q, name='scaled_free_counts')
    count_constraint = scaled_affine[:, 0] + scaled_affine[:, 1:] @ t >= 0
    primal_constraints = [count_constraint]
    for M in scaled_matrices:
        primal_constraints.append(M[0] + sum(M[j+1] * t[j] for j in range(q)) >> 0)
    primal = cp.Problem(cp.Minimize(0), primal_constraints)
    print('FRESH PRIMAL SOLVE', data['parameters'], solver, flush=True)
    primal.solve(solver=solver, verbose=True, **options)
    primal_record = {'status': primal.status, 'value': float(primal.value) if np.isfinite(primal.value) else None,
                     'statistics': statistics(primal)}
    if t.value is not None:
        primal_record['free_counts'] = (variable_scale * t.value).tolist()
        primal_record['scaled_count_minimum'] = float(np.min(scaled_affine @ np.r_[1., t.value]))
        primal_record['scaled_block_minimum_eigenvalues'] = [
            float(np.linalg.eigvalsh(np.einsum('j,jab->ab', np.r_[1., t.value], M))[0])
            for M in scaled_matrices]

    z = cp.Variable(len(affine), nonneg=True, name='scaled_count_multipliers')
    Ys = [cp.Variable(M.shape[1:], symmetric=True, name=f'Y_{b["H"]}_{b["ell"]}')
          for M, b in zip(scaled_matrices, data['blocks'])]
    contraction = [scaled_affine[:, j] @ z + sum(cp.sum(cp.multiply(Y, M[j]))
                   for Y, M in zip(Ys, scaled_matrices)) for j in range(q+1)]
    dual_constraints = [Y >> 0 for Y in Ys]
    dual_constraints += [contraction[j] == 0 for j in range(1, q+1)]
    dual_constraints += [cp.sum(z) + sum(cp.trace(Y) for Y in Ys) == 1]
    dual = cp.Problem(cp.Maximize(-contraction[0]), dual_constraints)
    print('FRESH DUAL SOLVE', solver, flush=True)
    dual.solve(solver=solver, verbose=True, **options)
    dual_record = {'status': dual.status, 'value': float(dual.value),
                   'statistics': statistics(dual)}
    directions = []
    if z.value is not None and all(Y.value is not None for Y in Ys):
        zv = np.array(z.value)
        Yv = [(Y.value + Y.value.T)/2 for Y in Ys]
        combined = scaled_affine.T @ zv + sum(
            np.einsum('ab,jab->j', Y, M) for Y, M in zip(Yv, scaled_matrices))
        dual_record.update({
            'scaled_count_multipliers': zv.tolist(),
            'original_count_multipliers': (zv / row_scales).tolist(),
            'combined_coefficients_in_1_scaled_xi': combined.tolist(),
            'max_coordinate_cancellation_residual': float(np.max(np.abs(combined[1:]))),
            'normalization_residual': float(abs(np.sum(zv) + sum(np.trace(Y) for Y in Yv)-1)),
            'minimum_count_multiplier': float(np.min(zv)),
            'blocks': []})
        for b, Y, scale in zip(data['blocks'], Yv, block_scales):
            vals, U = np.linalg.eigh(Y)
            dual_record['blocks'].append({'H': b['H'], 'ell': b['ell'],
                'scaled_Y': Y.tolist(), 'original_Y': (Y/scale).tolist(),
                'eigenvalues': vals.tolist(), 'block_scale': scale})
            for i, eta in enumerate(vals):
                if eta > 1e-10:
                    directions.append({'H': b['H'], 'ell': b['ell'],
                        'eigenvalue': float(eta), 'original_eigenvalue': float(eta/scale),
                        'f': U[:, i].tolist()})
        assert dual.value > 1e-8, 'No positive numerical dual objective found.'
        assert np.max(np.abs(combined[1:])) < 1e-7, 'Dual cancellation is not accurate enough.'
    else:
        raise RuntimeError('No dual solution returned')
    report = {'parameters': data['parameters'], 'model_sha256': hashlib.sha256(data_bytes).hexdigest(),
        'input': str(data_path.resolve()), 'input_provenance': 'Fresh exact article model only; no prior directions read.',
        'count_convention': data['count_convention'], 'profile_order': data['profile_order'],
        'variable_scaling': {'formula': 'xi_j = variable_scale * t_j', 'variable_scale': variable_scale},
        'count_row_scales': row_scales.tolist(), 'block_scales': block_scales,
        'scaling_convention': 'A_scaled = A * diag(1,s,...,s) / row_scale; M_scaled_j = M_j * coordinate_scale_j / block_scale. Original z=z_scaled/row_scale, original Y=Y_scaled/block_scale. No profile congruence.',
        'solver_options': options, 'primal': primal_record, 'dual': dual_record,
        'directions': directions, 'elapsed_seconds': time.time()-started,
        'proof_status': 'NUMERICAL ONLY; rational directions and exact nonnegative LP still required.'}
    output_path.write_text(json.dumps(report, indent=2) + '\n')
    print('RESULT', json.dumps({'primal_status': primal.status, 'dual_status': dual.status,
        'dual_objective': dual.value, 'direction_count': len(directions),
        'max_cancel_residual': dual_record['max_coordinate_cancellation_residual'],
        'elapsed_seconds': report['elapsed_seconds']}), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data', type=Path, default=HERE/'data.json')
    parser.add_argument('--output', type=Path, default=HERE/'numerical_solution.json')
    parser.add_argument('--solver', choices=('CLARABEL', 'SCS'), default='CLARABEL')
    parser.add_argument('--variable-scale', type=float, default=1e10)
    args = parser.parse_args()
    run(args.data, args.output, args.solver, args.variable_scale)
