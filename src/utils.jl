# Frequently-used static types
Mat3{T} = SMatrix{3, 3, T, 9} where T
Vec3{T} = SVector{3, T} where T

function determine_grid_size(lattice::AbstractMatrix, Ecut; kpoints=[[0, 0, 0]],
                             supersampling=2)
    # Lattice  and reciprocal lattice
    recip_lattice = SMatrix{3, 3}(lattice)
    recip_lattice = 2π * inv(Matrix(lattice'))
    kcart = [recip_lattice * k for k in kpoints]

    cutoff_qsq = 2 * supersampling^2 * Ecut
    # For a particular k-Point (in cartesian coordinates), the integer coordinates
    # [m n o] of the complementary reciprocal lattice vectors B satisfy
    #     |B * [m n o] + k|^2 ≤ cutoff_qsq
    # Now
    #     |B * [m n o] + k| ≥ abs(|B * [m n o]| - |k|) = |B * [m n o]| - |k|
    # provided that |k| ≤ |B|, which is typically the case. Therefore
    #     |[m n o]| / |B^{-1}| ≤ |B * [m n o]| ≤ sqrt(cutoff_qsq) + |k|
    # (where |B^{-1}| is the operator norm of the inverse of B), such that
    #     |[m n o]| ≤ (sqrt(cutoff_qsq) + |k|) * |B^{-1}|
    # In the extremal case, m = o = 0, such that
    #    n_max_trial = (sqrt(cutoff_qsq) + |k|) * |B^{-1}|
    #                = (sqrt(cutoff_qsq) + |k|) * |A| / 2π

    # Estimate trial upper bound n_max
    max_k = maximum(norm.(kcart))
    @assert max_k ≤ opnorm(recip_lattice)
    trial_n_max = ceil(Int, (max_k + sqrt(cutoff_qsq)) * opnorm(lattice) / 2π)
    # Determine actual n_max (trial_n_max is extended by 1 for safety)
    trial_n_range = -trial_n_max-1:trial_n_max+1
    n_max = 0
    for coord in CartesianIndices((trial_n_range, trial_n_range, trial_n_range))
        # TODO There certainly are more clever ways to do this than a triple loop ...
        energy(q) = sum(abs2, recip_lattice * q) / 2
        if any(energy([coord.I...] + k) ≤ supersampling^2 * Ecut for k in kcart)
            @assert all(abs.([coord.I...]) .<= trial_n_max)
            n_max = max(n_max, maximum(abs.([coord.I...])))
        end
    end
    return 2 * n_max + 1
end
