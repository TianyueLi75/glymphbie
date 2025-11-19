#ifndef _STOKESBIO_HPP_
#define _STOKESBIO_HPP_

#include "csbq.hpp"
// #include <slender_element.hpp>

// StokesBIE by Dhairya Malhotra, 2025: 
// https://github.com/dmalhotra/stokes-periodize

/**
 * PVFMM cannot handle combined field kernel. Compute SL and DL separately and add them.
 */
template <class Real> class StokesBIO {
  public:

    StokesBIO(const Real SL_scal, const Real DL_scal, const sctl::Comm comm);

    /**
     * Set periodicity.
     *
     * @param[in] periodicity periodicity type (NONE, X, XY, XYZ).
     *
     * @param[in] period_length length of the periodic box in each dimension.
     * Must be positive if periodicity is not NONE.
     *
     * @remark Periodicity only supported in 3D and with PVFMM.
     */
    void SetPeriodicity(sctl::Periodicity periodicity, Real period_length = 0);

    /**
     * Specify quadrature accuracy tolerance.
     *
     * @param[in] tol quadrature accuracy.
     */
    void SetAccuracy(Real tol);

    /**
     * Add an element-list.
     *
     * @param[in] elem_lst an object (of type ElemLstType, derived from the
     * base class ElementListBase) that contains the description of a list of
     * elements.
     *
     * @param[in] name a string name for this element list: note that this will be ordered alphabetically or numerically.
     */
    void AddElemList(const sctl::SlenderElemList<Real>& elem_lst, const std::string& name);

    /**
     * Delete an element-list.
     *
     * @param[in] name name of the element-list to return.
     */
    void DeleteElemList(const std::string& name);

    /**
     * Delete an element-list.
     */
    void DeleteElemList();

    /**
     * Set target point coordinates.
     *
     * @param[in] Xtrg the coordinates of target points in array-of-struct
     * order: {x_1, y_1, z_1, x_2, ..., x_n, y_n, z_n}
     */
    void SetTargetCoord(const sctl::Vector<Real>& Xtrg);

    /**
     * Set target point normals.
     *
     * @param[in] Xn_trg the coordinates of target points in array-of-struct
     * order: {nx_1, ny_1, nz_1, nx_2, ..., nx_n, ny_n, nz_n}
     */
    void SetTargetNormal(const sctl::Vector<Real>& Xn_trg);

    /**
     * Get local dimension of the boundary integral operator. Dim(0) is the
     * input dimension and Dim(1) is the output dimension.
     */
    sctl::Long Dim(sctl::Integer k) const;

    /**
     * Setup the boundary integral operator.
     */
    void Setup() const;

    /**
     * Clear setup data.
     */
    void ClearSetup() const;

    /**
     * Evaluate the boundary integral operator.
     *
     * @param[out] U the potential computed at each target point in
     * array-of-struct order.
     *
     * @param[in] F the charge density at each surface discretization node in
     * array-of-struct order.
     */
    void ComputePotential(sctl::Vector<Real>& U, const sctl::Vector<Real>& F) const;

    /**
     * Scale input vector by sqrt of the area of the element.
     * TODO: replace by sqrt of surface quadrature weights (not sure if it makes a difference though)
     */
    void SqrtScaling(sctl::Vector<Real>& U) const;

    /**
     * Scale input vector by inv-sqrt of the area of the element.
     * TODO: replace by inv-sqrt of surface quadrature weights (not sure if it makes a difference though)
     */
    void InvSqrtScaling(sctl::Vector<Real>& U) const;


  private:

    const sctl::Stokes3D_FxU ker_FxU;
    const sctl::Stokes3D_DxU ker_DxU;
    const sctl::Stokes3D_FxUP ker_FxUP;
    const sctl::Stokes3D_FSxU ker_FSxU;

    const sctl::Comm comm_;
    const Real SL_scal_, DL_scal_;
    sctl::BoundaryIntegralOp<Real, sctl::Stokes3D_FxU> LayerPotenSL;
    sctl::BoundaryIntegralOp<Real, sctl::Stokes3D_DxU> LayerPotenDL;
};

#endif