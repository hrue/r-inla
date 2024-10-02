{
	// this compute c0,c1,c2 at (eta,kappa,xi,alpha,y)

	double t1, t10, t102, t104, t108, t109, t11, t12, t123, t128, t13, t132, t14, t149, t15, t156, t158, t16,
		t161, t163, t164, t171, t173, t19, t20, t207, t21, t212, t22, t23, t25, t29, t3, t31, t32, t34,
		t35, t37, t38, t4, t42, t43, t44, t47, t48, t5, t50, t53, t54, t57, t6, t66, t69, t7, t73, t75,
		t76, t77, t8, t80, t82, t89, t9;
	double eeta = exp(eta);

	t1 = log(kappa);
        t3 = eeta;
        t7 = pow(alpha, 0.1e1 / kappa);
        t9 = pow(0.1e1 - t7, -xi);
        t10 = t9 - 0.1e1;
        t12 = 0.1e1 + y / t3 * t10;
        t13 = 0.1e1 / xi;
        t14 = pow(t12, -t13);
        t16 = log1p(- t14);
        t21 = log(xi * t3 / t10);
        t23 = log(t12);
        t25 = t1 + (kappa - 0.1e1) * t16 - t21 + t23 * (-t13 - 0.1e1);
	c0 = t25;

	t1 = y * kappa;
        t3 = pow(alpha, 0.1e1 / kappa);
        t4 = 0.1e1 - t3;
        t5 = pow(t4, xi);
        t7 = 1/eeta;
        t8 = eeta;
        t9 = pow(t4, -xi);
        t13 = 0.1e1 / xi;
        t14 = pow(t7 * (t8 + y * t9 - y), t13);
        t16 = t8 * t5;
        t20 = t14 * y;
        t31 = -(t1 - t1 * t5 + xi * t14 * t16 - xi * t8 * t5 - t20 + t20 * t5) * t13 / (t14 - 0.1e1) / (t16 + y - y * t5);
	c1 = t31;

        t1 = SQR(xi); 
        t3 = 1/eeta;
        t4 = eeta;
        t5 = 0.1e1 / kappa;
        t6 = pow(alpha, t5);
        t7 = 0.1e1 - t6;
        t8 = pow(t7, -xi);
        t11 = t3 * (t4 + y * t8 - y);
        t12 = 0.1e1 / xi;
        t13 = pow(t11, t12);
        t14 = t13 - 0.1e1;
        t15 = t14 * t14;
        t19 = pow(t7, xi);
        t20 = t4 * t19;
        t22 = t20 + y - y * t19;
        t23 = t22 * t22;
        t29 = pow(alpha, 0.2e1 * t5);
        t32 = pow(alpha, 0.3e1 * t5);
        t34 = pow(0.1e1 - 0.3e1 * t6 + 0.3e1 * t29 - t32, xi);
        t35 = t34 * xi;
        t37 = exp(0.2e1 * eta);
        t38 = kappa * t37;
        t42 = pow(0.1e1 - 0.2e1 * t6 + t29, xi);
        t43 = t13 * t42;
        t44 = xi * t37;
        t47 = pow(t11, 0.3e1 * t12);
        t48 = t47 * t42;
        t50 = t1 * t37;
        t53 = pow(t11, 0.2e1 * t12);
        t54 = t53 * t34;
        t57 = t53 * t42;
        t66 = t13 * t34;
        t69 = t47 * t34;
        t73 = t42 * xi;
        t75 = t53 * t19;
        t76 = SQR(y); 
        t77 = t76 * kappa;
        t80 = y * t4;
        t82 = t1 * t19;
        t89 = t35 * t38 + t43 * t44 + t48 * t44 + t48 * t50 + 0.3e1 * t54 * t50 - 0.2e1 * t57 * t44 + 0.2e1 * t54 * t44 - 0.3e1 * t57 * t50 + 0.3e1 * t43 * t50 - 0.3e1 * t66 * t50 - t69 * t44 - t69 * t50 - t66 * t44 - t73 * t38 - 0.3e1 * t75 * t77 - t75 * t80 - t82 * t80 + 0.3e1 * t77 * t13 * t19 + t20 * y * t13;
        t102 = t34 * t1;
        t104 = t42 * t1;
        t108 = xi * kappa;
        t109 = t108 * t4;
        t123 = t13 * t4 * t19;
        t128 = t53 * t76;
        t132 = t66 * t80 + t66 * t77 + 0.2e1 * t57 * t80 + 0.3e1 * t57 * t77 - 0.2e1 * t43 * t80 - 0.3e1 * t43 * t77 - t54 * t77 - t54 * t80 - t102 * t80 + 0.2e1 * t104 * t80 + 0.2e1 * t66 * y * t109 - t54 * y * t109 + 0.2e1 * t57 * y * t109 - 0.4e1 * t43 * y * t109 + 0.2e1 * y * xi * kappa * t123 - t75 * y * t109 - t128 + t76 * t13 + t102 * t37 - t104 * t37;
        t149 = t108 * t37;
        t156 = t47 * t19;
        t158 = xi * t4 * y;
        t161 = t1 * t4 * y;
        t163 = y * kappa;
        t164 = t163 * t4;
        t171 = t80 * t13;
        t173 = t128 * kappa + 0.3e1 * t75 * t76 - 0.3e1 * t19 * t76 * t13 - t77 * t13 + 0.3e1 * t43 * t76 - 0.3e1 * t57 * t76 + t54 * t76 - t66 * t76 - t108 * t20 * y + t54 * t149 + 0.2e1 * t43 * t149 - 0.2e1 * t66 * t149 - t57 * t149 + t156 * t158 + t156 * t161 + t75 * t164 - 0.2e1 * t75 * t158 - 0.3e1 * t75 * t161 + xi * t19 * t171;
        t207 = 0.3e1 * t82 * t171 - t66 * t164 + t69 * t161 + t69 * t158 + t66 * t158 - 0.2e1 * t54 * t158 - 0.3e1 * t54 * t161 + t54 * t164 + 0.4e1 * t57 * t158 + 0.6e1 * t57 * t161 - 0.2e1 * t57 * t164 + 0.2e1 * t43 * t164 - 0.2e1 * t48 * t158 - 0.2e1 * t48 * t161 - 0.2e1 * t43 * t158 - 0.6e1 * t43 * t161 + 0.3e1 * t66 * t161 - t35 * t164 + 0.2e1 * t73 * t164 - t163 * t123;
        t212 = -0.1e1 / t1 / t15 / t14 / t23 / t22 * y * (t89 + t132 + t173 + t207);
	c2 = t212;
}
