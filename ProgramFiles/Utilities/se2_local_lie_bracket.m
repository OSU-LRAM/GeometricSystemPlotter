function bracket = se2_local_lie_bracket(u,v)
% Calculate the lie bracket of two vectors in se(2) ( TeSE(2) )

	bracket = { (v{3}.*u{2}) - (u{3}.*v{2});
				(u{3}.*v{1}) - (v{3}.*u{1});
				zeros(size(u{3}))};

end