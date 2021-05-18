function[ dilated_points ] = dilatation( maillage , VertConn , order )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fonction permettant la dilatation d'une r�gion d'inter�t � l'ordre
% souhait�

% ENTREE: feature_points : les points du maillage que l'on souhaite dilater
%                           (marqu�s par un �tat haut)
%         VertConn : le vecteur de connectivit�
%         order : le nombre de dilattion successives que le porogramme va
%         effectuer

% SORTIE : dilated_points : un vecteur binaire donnant les point du
%                           maillage dilat�s

% REF: C. Rossl, L. Kobbelt, HP. Seidel, "Extraction of feature lines on
% triangulated surfaces using morphological operators"

% Date de cr�ation: vendredi 13 janvier 2006

%Script: Benoit Cottereau
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dilated_points = zeros( 1 , length( maillage ) );
feature_points = find( maillage );
acc = [];
for indice = 1 : length( feature_points )
    dilated_points( feature_points( indice ) ) = 1;
    dilated_points( VertConn{ feature_points( indice ) } ) = 1; % ainsi que son voisinage du premier ordre
end


if order > 1
   dilated_points = dilatation( dilated_points , VertConn , order - 1 ); 
end

return

