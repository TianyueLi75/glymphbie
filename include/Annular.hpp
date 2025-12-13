#ifndef _ANNULAR_HPP_
#define _ANNULAR_HPP_

// #include "csbq.hpp"
#include "csbq/slender_element.hpp"
#include "csbq/slender_element.cpp"


/**
 * @brief Implements a class for the annular space between two channels 
 *        using the Convergent slender-body quadrature (CSBQ): csbq.readthedocs.io
 *        Both serial and with MPI support using SCTL (sctl.readthedocs.io/en/latest/) 
 * 
 * @tparam Real: double or float.
 */
template <class Real> class Annular {

    public:
        /**
         * @brief Constructor (serial)
         * 
         * @param Xc_inner_ : 3D coordinates of the centerline of the inner channel, ordered [x1;y1;z1;x2;y2;z2;...] (for periodic channel, whole channel must be within unit box [0,1]^3)
         * @param Xc_outer_ : 3D coordinates of the centerline of the outer channel, same as above.
         * @param r_inner_ : vector of radius values of the inner channel, at each centerline node
         * @param r_outer_ : vector of radius values of the outer channel, at each centerline node
         */
        Annular(const sctl::Vector<Real> Xc_inner_, const sctl::Vector<Real> Xc_outer_, const sctl::Vector<Real> r_inner_, const sctl::Vector<Real> r_outer_);
        
        /**
         * @brief Constructor (with MPI). 
         * 
         * @param Xc_inner_ : 3D coordinates of the centerline of the inner channel, ordered [x1;y1;z1;x2;y2;z2;...] (for periodic channel, whole channel must be within unit box [0,1]^3)
         * @param Xc_outer_ : 3D coordinates of the centerline of the outer channel, same as above.
         * @param r_inner_ : vector of radius values of the inner channel, at each centerline node
         * @param r_outer_ : vector of radius values of the outer channel, at each centerline node
         * @param comm_ : MPI Communicator in the SCTL library
         */
        Annular(const sctl::Vector<Real> Xc_inner_, const sctl::Vector<Real> Xc_outer_, const sctl::Vector<Real> r_inner_, const sctl::Vector<Real> r_outer_, sctl::Comm comm_);
        
        /**
         * @brief Destructor.
         * 
         */
        ~Annular();

        /**
         * @brief Set the centerline nodes for the inner channel
         * 
         * @param X_in 
         */
        void SetInnerXc(const sctl::Vector<Real> X_in);

        /**
         * @brief Set the centerline nodes for the outer channel
         * 
         * @param X_in 
         */
        void SetOuterXc(const sctl::Vector<Real> X_in);

        /**
         * @brief Set the radii for the inner channel
         * 
         * @param r_in 
         */
        void SetInnerR(const sctl::Vector<Real> r_in);

        /**
         * @brief Set the radii for the outer channel
         * 
         * @param r_in 
         */
        void SetOuterR(const sctl::Vector<Real> r_in);

        /**
         * @brief Set up the annular geometry. 
         * 
         * Using centerline and radius information, creates circular cross-sectional channels with the number
         * of Fourier nodes specified.
         * 
         * @param Nelem_ : Total number of panels along axial direction
         * @param ElemOrder_ : Order of Chebyshev quadrature nodes on each panel
         * @param FourierOrder_ : Number of Fourier nodes in azimuthal direction
         * @param drdt_inner : time derivatve of radius of the inner wall, for slip velocity calculations
         * @param drdt_outer : time derivatve of radius of the outer wall, for slip velocity calculations
         * @param force_setup : If true, recompute the centerline nodes of specified order using interpolation 
         *                      even if the setup booleans are true.
         * @return std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> 
         */
        std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> Setup(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_inner, sctl::Vector<Real> drdt_outer, bool force_setup = false);
        
        /**
         * @brief Set up the annular geometry, with MPI support. 
         * 
         * Using centerline and radius information, creates circular cross-sectional channels with the number
         * of Fourier nodes specified. Divide the total number of panels evenly among the available MPI processes. 
         * Only the local segment information is stored.
         * 
         * @param Nelem_ : Total number of panels along axial direction
         * @param ElemOrder_ : Order of Chebyshev quadrature nodes on each panel
         * @param FourierOrder_ : Number of Fourier nodes in azimuthal direction
         * @param drdt_inner : time derivatve of radius of the inner wall, for slip velocity calculations
         * @param drdt_outer : time derivatve of radius of the outer wall, for slip velocity calculations
         * @param force_setup : If true, recompute the centerline nodes of specified order using interpolation 
         *                      even if the setup booleans are true.
         * @return std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> 
         */
        std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> Setup_mpi(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_inner_, sctl::Vector<Real> drdt_outer_, bool force_setup = false);
        
        /**
         * @brief Set up the inner channel
         * 
         * If only the inner channel is changed, with different slip conditions, or have a different order, 
         * only SetupInner should to be called.
         * 
         * @param Nelem_ : Total number of panels along axial direction
         * @param ElemOrder_ : Order of Chebyshev quadrature nodes on each panel
         * @param FourierOrder_ : Number of Fourier nodes in azimuthal direction
         * @param drdt_inner : time derivatve of radius of the inner wall, for slip velocity calculations
         * @param force_setup : If true, recompute the centerline nodes of specified order using interpolation 
         *                      even if the setup booleans are true.
         * @return std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> 
         */
        std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> SetupInner(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_inner, bool force_setup = false);
        
        /**
         * @brief Set up the outer channel
         * 
         * If only the outer channel is changed, with different slip conditions, or have a different order, 
         * only SetupOuter should to be called.
         * 
         * @param Nelem_ : Total number of panels along axial direction
         * @param ElemOrder_ : Order of Chebyshev quadrature nodes on each panel
         * @param FourierOrder_ : Number of Fourier nodes in azimuthal direction
         * @param drdt_outer : time derivatve of radius of the inner wall, for slip velocity calculations
         * @param force_setup : If true, recompute the centerline nodes of specified order using interpolation 
         *                      even if the setup booleans are true.
         * @return std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> 
         */
        std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> SetupOuter(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_outer, bool force_setup = false);
        
        std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> SetupInner_mpi(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_inner, bool force_setup = false);
        
        std::tuple<sctl::Vector<Real>, sctl::Vector<Real>, sctl::Vector<Real>> SetupOuter_mpi(const sctl::Long Nelem_, const sctl::Long ElemOrder_, const sctl::Long FourierOrder_, sctl::Vector<Real> drdt_outer, bool force_setup = false);
        
        /**
         * @brief Get the inner channel node coordinates and normal vectors
         * 
         * @param X_out : location to load node coordinates, order [x1;y1;z1;x2;y2;z2;...]
         * @param Xn_out : location to load normal vectors.
         */
        void GetInnerCoord(sctl::Vector<Real>* X_out, sctl::Vector<Real>* Xn_out);

        void GetOuterCoord(sctl::Vector<Real>* X_out, sctl::Vector<Real>* Xn_out);

        /**
         * @brief Get all channel node coordinates and normal vectors, order: X_inner then X_outer
         * 
         * @param X : location to load node coordinates, order [x1;y1;z1;x2;y2;z2;...]
         * @param Xn : location to load normal vectors.
         */
        void GetNodeCoord(sctl::Vector<Real>* X, sctl::Vector<Real>* Xn); 

        /**
         * @brief Get the 3D centerline coordinates on a straight line at (y=0.5,z=0.5)
         *        
         *        Panels are set equidistant with size 1/<Nelem_>. On each, x nodes are order <ElemOrder_> Chebyshev nodes
         *        More complicated centerlines can be achieved by taking only the x coordinates and resetting y and z.
         * 
         * @param Nelem_ : Total number of panels
         * @param ElemOrder_ : Chebyshev order on each panel.
         * @return sctl::Vector<Real> 
         */
        static sctl::Vector<Real> GetCenterLine(const sctl::Long Nelem_, const sctl::Long ElemOrder_);

        /**
         * @brief Get the 3D centerline coordinates on a straight line at (y=0.5,z=0.5), MPI support
         * 
         *        Split the total number of panels <Nelem_> evenly by the number of MPI processes first,
         *        only the local segment of the centerline is returned. 
         *        The number of panels on current process and the global indexing of the first panel on this process are returned as well.
         * 
         * @param Nelem_ 
         * @param ElemOrder_ 
         * @param comm_ 
         * @return std::tuple<sctl::Long, sctl::Long, sctl::Vector<Real>> 
         */
        static std::tuple<sctl::Long, sctl::Long, sctl::Vector<Real>> GetCenterLine_mpi(const sctl::Long Nelem_, const sctl::Long ElemOrder_, sctl::Comm comm_);
        
        /**
         * @brief Returns the minimum radius *over local segment if using MPI* of the inner channel
         *
         *        This value can be used for scaled (convergent) BIO for slender-body objects.
         * 
         * @return Real 
         */
        Real GetMinRadius();

        /**
         * @brief Determines if a given point is in between the *local segments if using MPI* channels
         * 
         * @param X_ : Vector of target to determine
         * @return sctl::Vector<sctl::Long> 
         */
        sctl::Vector<sctl::Long> InDomain(sctl::Vector<Real> X_); // given a vector of x positions, return which nodes are in bewteen outer and inner channels V.S. not.
        
        sctl::Vector<sctl::Long> InDomain_mpi(sctl::Vector<Real> X_);

        /**
         * @brief Get the Inner SlenderElemList object for the inner channel
         * 
         * @see SlenderElemList
         * 
         * @return sctl::SlenderElemList<Real> 
         */
        sctl::SlenderElemList<Real> GetInnerElemList();

        sctl::SlenderElemList<Real> GetOuterElemList();

        /**
         * @brief Evaluate the time derivative on surface nodes of both the inner and outer channel radii from the given dr/dt at centerline nodes.
         * 
         * @return sctl::Vector<Real> 
         */
        sctl::Vector<Real> GetVslip();

        sctl::Vector<Real> GetVslip_mpi();

        sctl::Vector<Real> GetVslipInner();

        sctl::Vector<Real> GetVslipOuter();

        sctl::Vector<Real> GetVslipInner_mpi();

        sctl::Vector<Real> GetVslipOuter_mpi();

        /**
         * @brief Writes the inner and outer channel geometry along with the given fields <F_inner>, <F_outer>, into vtk format, using the WriteVTK functions in SlenderElemList.
         *
         * @see SlenderElemList
         * 
         * @param filename_inner 
         * @param filename_outer 
         * @param F_inner : Field to be written on the inner channel element_list nodes
         * @param F_outer : Field to be written on the outer channel element_list nodes
         * @param comm_ : SCTL communicator
         */
        void WriteVTK(std::string filename_inner, std::string filename_outer, sctl::Vector<Real> F_inner, sctl::Vector<Real> F_outer, sctl::Comm comm_);

        // TODO: setup adaptive grid, and return new centerline after each solve.



    private:
        sctl::Vector<Real> _Xc_inner, _Xc_outer, _r_inner, _r_outer; // input values
        sctl::Vector<Real> _X_inner, _X_outer, _Xn_inner, _Xn_outer; // computed values from CSBQ; NOTE: Normal points into fluid.
        Real _min_radius = 0.; // computed value
        sctl::Long _Nelem_inner, _Nelem_outer, _ElemOrder_inner, _ElemOrder_outer, _FourierOrder_inner, _FourierOrder_outer; 
        sctl::Vector<Real> _drdt_inner, _drdt_outer;

        // MPI specific params
        sctl::Comm _comm;
        sctl::Long _loc_elem_dsp_inner, _loc_elem_dsp_outer;
        sctl::Long _loc_elem_cnt_inner, _loc_elem_cnt_outer;

        sctl::SlenderElemList<Real> _elem_lst_outer, _elem_lst_inner;

        bool _SetupInner_bool = false;
        bool _SetupOuter_bool = false;

        void InterpR(sctl::Vector<Real>& trg_r, const sctl::Vector<Real> src_r, const sctl::Vector<Real> src_x, const sctl::Vector<Real> trg_x);
};









#endif