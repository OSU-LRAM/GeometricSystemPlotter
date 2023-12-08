%Plots decides which plotting function to use based on salp definition
function plotSalpFrame(ax,salp,pose,shape,forcedef,t)

    %If force is not defined, don't try to pass an undefined variable
    if nargin < 5

        switch salp.geometry.type
            %If salp is a series of rigid links plot it that way
            case 'n-link chain'
                plotSalpFrame_NLinkChain(ax,salp,pose,shape);
            %If salp is continuous plot it that way
            case 'general curvature'
                plotSalpFrame_GeneralCurvature(ax,salp,pose,shape);
            otherwise
                error('Incorrectly defined salp geometry type');
        end

    else

        switch salp.geometry.type
            case 'n-link chain'
                plotSalpFrame_NLinkChain(ax,salp,pose,shape,forcedef,t);
            case 'general curvature'
                plotSalpFrame_GeneralCurvature(ax,salp,pose,shape,forcedef,t);
            otherwise
                error('Incorrectly defined salp geometry type');
        end

    end

end