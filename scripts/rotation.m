classdef rotation
    %ROTATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        %Property1
    end
    
    methods (Static)
        function obj = rotation()%inputArg1,inputArg2)
            %ROTATION Construct an instance of this class
            %   Detailed explanation goes here
            %obj.Property1 = inputArg1 + inputArg2;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [uCross] = crossProductEquivalent(u)
            % crossProductEquivalent : Outputs the cross-product-equivalent
            %                          matrix uCross such that for arbitrary
            %                          3-by-1 vectors u and v:
            %                          cross(u,v) = uCross*v.
            %
            % INPUTS
            %
            % u ---------- 3-by-1 vector
            %
            %
            % OUTPUTS
            %
            % uCross ----- 3-by-3 skew-symmetric cross-product equivalent matrix
            
            uCross = [     0, -u(3),  u(2) ;
                        u(3),     0, -u(1) ;
                       -u(2),  u(1),    0 ];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [quatOut] = quatmult_S(quat1,quat2)
            %QUATMULT_S Computes the Schuster quaternion product of two
            %quaternions
            %   Computes qOut = q1 * q2, where * is the quaternion product
            %   defined by Malcolm Schuster, which is the composition formula
            %   for rotating-space successive rotations.
            %
            % Reference: ROTATIONS, TRANSFORMATIONS, LEFT QUATERNIONS, RIGHT
            %            QUATERNIONS? - Renato Zanetti
            
            % TODO: double check order of q1, q2
            
            if length(quat2) ~= 4,
                disp("incorrect input");
            end
            if length(quat1) ~= 4,
                disp("incorrect input");
            end
            
            if ~iscolumn(quat2),
                quat2 = transpose(quat2);
            end
            if ~iscolumn(quat1),
                quat1 = transpose(quat1);
            end
            
            %quat1 = rotation.normedQuat(quat1);
            %quat2 = rotation.normedQuat(quat2);
            
            q1w = quat1(1); % scalar
            q2w = quat2(1); % scalar
            
            q1vec = quat1(2:4); % vector
            q2vec = quat2(2:4); % vector
            
            qOut_w = q1w*q2w - dot(q1vec, q2vec);
            qOut_vec = q1w*q2vec + q2w*q1vec - cross(q1vec, q2vec);
            
            quatOut = [qOut_w; qOut_vec];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [quatOut] = quatmult_H(quat2,quat1)
            %QUATMULT_H Computes the Hamilton quaternion product of two
            %quaternions
            %   Computes qOut = q2 * q1, where * is the quaternion product
            %   defined by William Hamilton, which is the composition formula
            %   for fixed-space successive rotations.
            %
            % Reference: ROTATIONS, TRANSFORMATIONS, LEFT QUATERNIONS, RIGHT
            %            QUATERNIONS? - Renato Zanetti
            
            if length(quat2) ~= 4,
                disp("incorrect input");
            end
            if length(quat1) ~= 4,
                disp("incorrect input");
            end
            
            if ~iscolumn(quat2),
                quat2 = transpose(quat2);
            end
            if ~iscolumn(quat1),
                quat1 = transpose(quat1);
            end
            
            %quat1 = rotation.normedQuat(quat1);
            %quat2 = rotation.normedQuat(quat2);
            
            q1w = quat1(1); % scalar
            q2w = quat2(1); % scalar
            
            q1vec = quat1(2:4); % vector
            q2vec = quat2(2:4); % vector
            
            qOut_w = q1w*q2w - dot(q2vec, q1vec);
            qOut_vec = q1w*q2vec + q2w*q1vec + cross(q2vec, q1vec);
            
            quatOut = [qOut_w; qOut_vec];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [quatOut] = normedQuat(quat)
            %NORMEDQUAT Returns the normalised input quaternion
            %
            quatOut = quat ./ norm(quat);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Tmat] = quat2Tmat(quat)
            %QUAT2TMAT Returns the (passive) transformation matrix corresponding
            %to the input quaternion
            %
            qw = quat(1);
            qvec = quat(2:4);
            qcross = rotation.crossProductEquivalent(qvec);
            Tmat = eye(3) - 2*qw*qcross + 2*qcross^2;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Rmat] = quat2Rmat(quat)
            %QUAT2RMAT Returns the (active) rotation matrix corresponding to the
            %input quaternion
            %
            
            if ~iscolumn(quat),
                quat = transpose(quat);
            end
            
            qw = quat(1);
            qvec = quat(2:4);
            qcross = rotation.crossProductEquivalent(qvec);
            Rmat = eye(3) + 2*qw*qcross + 2*qcross^2;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [quat] = angleaxis2quat(theta, nHat)
            %ANGLEAXIS2QUAT Returns the quaternion associated with the (active)
            %axis-angle rotation of angle theta (in rad) about unit vector nHat
            %
            
            if ~iscolumn(nHat),
                nHat = transpose(nHat);
            end
            
            qw = cos(theta/2);
            qvec = sin(theta/2)*nHat;
            
            quat = [qw; qvec];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Tmat] = angleaxis2Tmat(theta, nHat)
            %ANGLEAXIS2TMAT Returns the (passive) transformation matrix
            %associated with the (active) axis-angle rotation of angle theta
            %(in rad) about unit vector nHat
            %   
            Tmat = rotation.quat2Tmat(rotation.angleaxis2quat(theta, nHat));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [Rmat] = angleaxis2Rmat(theta, nHat)
            %ANGLEAXIS2TMAT Returns the (active) rotation matrix
            %associated with the (active) axis-angle rotation of angle theta
            %(in rad) about unit vector nHat
            %   
            Rmat = rotation.quat2Rmat(rotation.angleaxis2quat(theta, nHat));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [e] = dcm2euler_312(R_BW)
            % dcm2euler : Converts a direction cosine matrix R_BW to Euler
            %             angles phi = e(1), theta = e(2), and psi = e(3)
            %             (in radians) for a 3-1-2 rotation.
            %
            % Let the world (W) and body (B) reference frames be initially
            % aligned. In a 3-1-2 order, rotate B away from W by angles psi
            % (yaw, about the body Z axis), phi (roll, about the body X axis),
            % and theta (pitch, about the body Y axis). R_BW can then be used to
            % cast a vector expressed in W coordinates as a vector in B
            % coordinates: vB = R_BW * vW
            %
            % INPUTS
            %
            % R_BW ------- 3-by-3 direction cosine matrix
            %
            %
            % OUTPUTS
            %
            % e ---------- 3-by-1 vector containing the Euler angles in radians:
            %              phi = e(1), theta = e(2), and psi = e(3)
            %
            %+-----------------------------------------------------------------+
            % References: D. Mellinger, N. Michael, and V. Kumar,
            %            "Trajectory generation and control for precise
            %            aggressive maneuvers with quadrotors,” The
            %            International Journal of Robotics Research,
            %            vol. 31, no. 5, pp. 664–674, 2012.
            
            phi = asin( R_BW(2,3) );
            
            theta = atan2( -R_BW(1,3), R_BW(3,3) );
            %theta = asin( -R_BW(1,3) / cos(phi) );
            
            psi = atan2( -R_BW(2,1), R_BW(2,2) );
            
            e = [phi; theta; psi];
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function [R_BW] = euler2dcm_312(e)
            % euler2dcm : Converts Euler angles phi = e(1), theta = e(2), and
            %             psi = e(3) (in radians) into a direction cosine matrix
            %             for a 3-1-2 rotation.
            %
            % Let the world (W) and body (B) reference frames be initially
            % aligned. In a 3-1-2 order, rotate B away from W by angles psi
            % (yaw, about the body Z axis), phi (roll, about the body X axis),
            % and theta (pitch, about the body Y axis). R_BW can then be used to
            % cast a vector expressed in W coordinates as a vector in B
            % coordinates: vB = R_BW * vW
            %
            % INPUTS
            %
            % e ---------- 3-by-1 vector containing the Euler angles in radians:
            %              phi = e(1), theta = e(2), and psi = e(3)
            %
            %
            % OUTPUTS
            %
            % R_BW ------- 3-by-3 direction cosine matrix
            
            phi   = e(1);
            theta = e(2);
            psi   = e(3);
            
            R_BW = rotation.angleaxis2Tmat(theta, [0 1 0]) * ...
                   rotation.angleaxis2Tmat(  phi, [1 0 0]) * ...
                   rotation.angleaxis2Tmat(  psi, [0 0 1]);
        end
        
    end
end

